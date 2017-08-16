#!/usr/bin/env python

###########
# Extractor for IMGT template sequences
########
import sys
import os.path
import time
import csv
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.support.ui import Select

import re

# selenium is necessary to deal with JS (<select> statements)

def delay_exec(statement):
    while True:
        try:
            exec(statement)
        except NoSuchElementException:
            #print("wait")
            time.sleep(0.1) # wait a bit :-)
        else:
            break
def delay_eval(statement):
    result = None
    while True:
        try:
            result = eval(statement)
        except NoSuchElementException:
            #print("wait")
            time.sleep(0.1) # wait a bit :-)
        else:
            break
    return(result.text) 

def set_option(driver, field, value):
    # searches the input element with name "field" and sets value "value"
    #print("setting field: " + field + " with value: " + value)
    element = Select(driver.find_element_by_name(field))
    #if element.is_multiple(): # not working
        #element.deselect_all()
    element.select_by_value(value)
    #select.selectByVisibleText("Value1");
    #element = driver.find_element_by_name(gdb_species)


# OPTIONS
#################
# select the folder where primers should be stored:
# select the species:
#############
# regular exressions and paths
# load 
print("call:")
print(sys.argv)
if len(sys.argv) != 7: # two args are required ([0] is name of script)
    sys.exit("Stopping because required args were not supplied: <phantomJS_binary> <out_loc_leader> <out_loc_exon> <species> <locus> <function>")
leader_file = sys.argv[2] # read from 2nd arg
exon_file = sys.argv[3]
if not os.path.exists(os.path.dirname(leader_file)):
    #print(os.path.dirname(leader_file))
    os.makedirs(os.path.dirname(leader_file))
# need to set capabilities for phantomjs -> change useragent to get results
dcap = dict(webdriver.DesiredCapabilities.PHANTOMJS)
dcap["phantomjs.page.settings.userAgent"] = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_8_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/29.0.1547.57 Safari/537.36"
phantom_executable = sys.argv[1] # read from 1st arg

driver = webdriver.PhantomJS(desired_capabilities = dcap, 
                            executable_path=phantom_executable,
                            service_log_path = os.path.devnull) #"/var/log/phantomjs/ghostdriver.log") # use phantomjs, store log in /var/
#driver = webdriver.Firefox() # opens browser for validation
driver.maximize_window()
base_url = "http://imgt.org/genedb/"
#########
#INPUT PARAMS: species (4), locus (5), function (6)
#########
# Identification

gdb_species = 'model.gene.id.species'
species = sys.argv[4]
gdb_molecular_component = "model.molComponent"
molecular_component = "any"
gdb_function = "model.allele.fcode"
function = sys.argv[6]
# Localization
gdb_locus = "model.locusLike"
locus = sys.argv[5]
# Classification
gdb_imgt_group = "model.groupLike" 
imgt_group = "any" # "any" is default
gdb_imgt_subgroup = "model.subgroup"
imgt_subgroup = "-1" # default (any)
input_options = {gdb_species : species,  gdb_function: function, gdb_locus: locus, gdb_imgt_group: imgt_group, gdb_imgt_subgroup: imgt_subgroup}
#####################
driver.get(base_url)
# set parameters:
for field in input_options:
    value = input_options[field]
    set_option(driver, field, value)
# molecular component: setting isn't necessary and would need to be hardcoded
#inputField = driver.find_element_by_name(gdb_molecular_component)
#inputField.send_keys(molecular_component)
# submit query
submit_button = driver.find_element_by_id("resultPage_0")
submit_button.submit()
# RESULTS PAGE
# sanity check: are there any results to retrieve?
try:
    driver.find_element_by_xpath("//*[contains(@class, 'noresult')]") # any element of class 'noresult' on page?
except NoSuchElementException: # this exception should be raised if results are available -> we can continue
    pass
else: # no exception occurred -> there weren't any results :-(
    sys.exit()

checkbox = driver.find_element_by_xpath("/html/body/div/form/p/input")
checkbox.click() # select all genes
# 1. extract leaders
checkbox = driver.find_element_by_xpath("/html/body/div/form/table[2]/tbody/tr[2]/td[3]/input[1]")
checkbox.click() 
set_option(driver, "selectedLabels", "L-PART1+L-PART2")
submit_button = driver.find_element_by_xpath("/html/body/div/form/input")
submit_button.click() # submit query
print("Retrieved leaders successfully :-)")
# store data
output_handle = open(leader_file, "w")
template_info = driver.find_element_by_xpath("/html/body/pre[2]").text # extract content from second <pre>
output_handle.write(template_info)
output_handle.close()
# 2. extract leader + variable region (L-PART1 + L-PART2 is extracted, although label says something different) + exon
driver.back()  # go back to previous browser page
# clear previous region selection and select new region
element = Select(driver.find_element_by_name("selectedLabels"))
element.deselect_by_value("L-PART1+L-PART2") 
set_option(driver, "selectedLabels", "L-PART1+V-EXON")
submit_button = driver.find_element_by_xpath("/html/body/div/form/input") # update submit button
submit_button.click() # submit query
print("Retrieved exons successfully :-)")
# store data
output_handle = open(exon_file, "w")
template_info = driver.find_element_by_xpath("/html/body/pre[2]").text # extract content from second <pre>
output_handle.write(template_info)
output_handle.close()
