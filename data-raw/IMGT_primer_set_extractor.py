###########
# IMGT PRIMER SET EXTRACTOR SCRIPT
############

import mechanize
import os.path
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
###################
# AMBIGUITY TABLE
####################
IUPAC_AMBIGUITIES = {"a":"a", "c" : "c", "g":"g", "t":"t", "m": {"a", "c"},\
                    "r" : {"a", "g"}, "w" : {"a", "t"}, "s" : {"c", "g"}, \
                    "y" : {"c", "t"}, "k" : {"g", "t"}, "v" : {"a", "c", "g"},\
                    "k" : {"g", "t"}, "v" : {"a", "c", "g"}, \
                    "h" : {"a", "c", "t"}, "d" : {"a", "g","t"}, \
                    "b" : {"c", "g", "t"},"n" : {"a", "c", "g", "t"}\
                    }
################
# OPTIONS
#################
# select the imgt groups for which we want to retrieve primer sets:
IMGT_groups = ["IGHV", "IGKV", "IGLV"]
# select the folder where primer fasta should be stored:
primer_folder = "IMGT_data/IMGT_primers"
# input for IMGT web interface:
# select the species:
species = "Homo sapiens"
species_url = species.replace(" ", "+")
# select the direction of the primers:
direction = "Sense"
#############
# regular exressions and paths
marker = "SetViewer_PrDB.pl"
flat_url = "http://www.imgt.org/IMGTPrimerDB/onePrimerToFlatFile.pl"
primer_marker = "Check_PrDB.pl" 
ID_expr = r"DE\t(.+)"
seq_expr = r"SQ\t.+\n\t+([a-zA-Z\(\)/\s]+)"
ambig_expr = r"\(([a-z/)]+)\)"
######
if not os.path.exists(primer_folder):
    os.makedirs(primer_folder)
for (group_idx,cur_group) in enumerate(IMGT_groups):
    #if cur_group != "IGLV":
        #continue
    print("Group: " + cur_group)
    ID_format_expr = "(" + cur_group + "\S*)"
    group_folder = os.path.join(primer_folder, cur_group)
    if not os.path.exists(group_folder):
        os.makedirs(group_folder)
    url = "http://www.imgt.org/IMGTPrimerDB/ListingPrimer_PrDB.pl?choice=set&SetGroupList=" + cur_group + "&SetEspeceList=" + species_url + "&SetDirection=" + direction + "&FormSet=Display+set+list"
    br = mechanize.Browser()
    br.set_handle_robots(False) # ignore robots
    page = br.open(url)
    links = br.links() # all links on the page
    sel_links = [l for l in links if marker in l.url]

    for (i,link) in enumerate(sel_links):
        #print link.text, link.url
        # follow the link :-)
        set_accession = link.text
        print("Set: " + set_accession)
        response = br.follow_link(link)
        html_location = "mechanize_results_" + str(i) + ".html"
        #with open(html_location, "w") as f:
            #f.write(response.read())
        #if link.text != "IPS000245":
        #    continue
        set_folder = os.path.join(group_folder, set_accession)
        #if not os.path.exists(set_folder):
            #os.makedirs(set_folder)
        # for every primer, retrieve the flat file
        set_links = br.links()
        sel_set_links = [l for l in set_links if primer_marker in l.url]
        set_html = response.read()
        #ref_marker = "http://www\.ncbi\.nlm\.nih\.gov:80"
        m = re.search(r".*?PMID.*?ExternalWin\('(" + ".*?)'\)", set_html, re.MULTILINE|re.DOTALL)
        b_ref = mechanize.Browser()
        b_ref.set_handle_robots(False) # ignore robots
        ref_result = b_ref.open(m.group(1))
        #meta name="author" content="Marks JD , et al." /><meta name="description" content="J Mol Biol. 1991 Dec 5;222(3):581-97. Research Support, Non-U.S. Gov't"
        ncbi_result = ref_result.read()
        #print(ncbi_result)
        ncbi_author = r"meta name=\"author\" content=\"([a-zA-Z]+)"
        m = re.search(ncbi_author, ncbi_result)
        NCBI_author = m.group(1)
        ncbi_year = r"meta name=\"description\" content=\".*?([0-9]+)"
        m = re.search(ncbi_year, ncbi_result)
        NCBI_year = m.group(1)
        fasta_name = set_accession + "_" + NCBI_author + str(NCBI_year)
        fasta_location = os.path.join(group_folder, fasta_name + ".fasta")
        set_ids = []
        set_seqs = []
        for (j, p_link) in enumerate(sel_set_links): # for every primer in cur set
            p_url = flat_url + "?Numacc=" + p_link.text # flat file url for cur primer
            p_loc = os.path.join(set_folder, "primer_" + str(j) + "_" + p_link.text + ".txt")
            primer_br = mechanize.Browser()
            primer_res = primer_br.open(p_url)
            primer_info = primer_res.read()
            m = re.search(r".*\<pre\>(.+)\</pre\>", primer_info, re.MULTILINE|re.DOTALL)
            if m:
                primer_info = m.group(1)
                m = re.search(ID_expr, primer_info)
                ID = m.group(1)
                # format ID (shorten):
                print(ID)
                m = re.search(ID_format_expr, ID)
                if m is None: # no detailed info available on the group
                    ID = cur_group
                else:
                    ID = m.group(1)
                ID = NCBI_author + str(NCBI_year) + "|" + ID + "|" + str(j+1) + "_fw"
                print(ID)
                set_ids.append(ID)
                m = re.search(seq_expr, primer_info)
                seq = m.group(1)
                seq = "".join(seq.split())
                # replace ambiguities of the form (a/b/c/d) in the sequence
                IUPAC_AMBIGUITIES.keys()[IUPAC_AMBIGUITIES.values().index("a")]
                print(seq)
                m = re.findall(ambig_expr, seq)
                iupac_nts = [IUPAC_AMBIGUITIES.keys()[IUPAC_AMBIGUITIES.values().index(set(amb.split("/")))] for amb in m]
                # substitute converted NTs in sequence
                for nt in iupac_nts:
                    seq = re.sub(ambig_expr, nt, seq, count=1)
                set_seqs.append(seq)
                #with open(p_loc, "w") as f:
                    #f.write(primer_info)

        # one set finished: write the fasta file
        output_handle = open(fasta_location, "w")
        #r = SeqRecord(Seq(set_seqs[0],IUPAC.ambiguous_dna))
        records = [SeqRecord(Seq(s,IUPAC.ambiguous_dna), \
                    id=set_ids[x], description="") for (x,s) in enumerate(set_seqs)]
        SeqIO.write(records, output_handle, "fasta")
        output_handle.close()



            

