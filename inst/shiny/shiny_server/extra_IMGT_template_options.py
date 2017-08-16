#!/usr/bin/env python

import os.path
import time
import csv
import sys
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

def get_option(driver, field):
    #print("getting field options " + field)
    element = Select(driver.find_element_by_name(field))
    options = []
    for opt in element.options:
        val = opt.text
        options.append(val)
    return(options)

# OPTIONS
#################
# select the folder where primers should be stored:
# select the species:
#############
# regular exressions and paths
# load 
if len(sys.argv) != 3: # two args are required ([0] is name of script)
    sys.exit("Stopping because required args were not supplied: <phantomJS_binary> <out_loc>")
options_folder = sys.argv[2] # read from 2nd arg
groups = ["IGHV", "IGKV", "IGLV"]

if not os.path.exists(options_folder):
    os.makedirs(options_folder)
# need to set capabilities for phantomjs -> change useragent to get results
dcap = dict(webdriver.DesiredCapabilities.PHANTOMJS)
dcap["phantomjs.page.settings.userAgent"] = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_8_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/29.0.1547.57 Safari/537.36"
#phantom_executable = '/home/mdoering/Downloads/phantomjs-2.1.1-linux-x86_64/bin/phantomjs'
phantom_executable = sys.argv[1] # read from 1st arg
driver = webdriver.PhantomJS(desired_capabilities = dcap, 
                            executable_path=phantom_executable,
                            service_log_path = os.path.devnull) # without browser window -> needs to be installed
print("Driver initialized successfully :-)")
#driver = webdriver.Firefox() # opens browser for validation
driver.maximize_window()
base_url = "http://imgt.org/genedb/"
#########
#INPUT PARAMS
#########
# Identification
gdb_species = 'model.gene.id.species'
gdb_molecular_component = "model.molComponent"
gdb_function = "model.allele.fcode"
# Localization
gdb_locus = "model.locusLike"
# Classification
gdb_imgt_group = "model.groupLike" 
gdb_imgt_subgroup = "model.subgroup"
input_options = [gdb_species, gdb_function, gdb_locus, gdb_imgt_group, gdb_imgt_subgroup]
# QUERY PAGE
driver.get(base_url)
# get parameters:
for field in input_options:
    options = get_option(driver, field)
    out_loc = options_folder + "/" + field + ".txt"
    f = open(out_loc, "w")
    options_str = '\n'.join(options) + '\n'
    f.writelines(options_str)



