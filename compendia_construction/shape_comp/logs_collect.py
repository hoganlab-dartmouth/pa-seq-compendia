#------------------------------------------
#
# Script Name: logs_collect.py
#
# Description: Collects run info from sra_to_salmon runs into a csv
#
# Args: 
#      dir 
#
# Author: Georgia Doing
#
# Date Created: 2021-02-11
#
# Email: Georgia.Doing.Gr@Dartmouth.edu
#
#-----------------------------------------
#
# Notes:
#
#
#
#-------------------------------------------

import sys
import glob
import csv
import json 

dirf = sys.argv[1]
log_files = glob.glob(dirf + "/*/*/*.salmon/aux_info/meta_info.json")
p_log_files = glob.glob(dirf + "/*/*/processes.log")
print(len(log_files))
print(len(p_log_files))
log_data = {}
for log in range(len(log_files)):
    #print(log)
    with open(log_files[log],'r') as lg:
        data = json.load(lg)
#        print(data.keys())
    ind = log_files[log][log_files[log].find('RP')-1:log_files[log].find('.salmon')]
    #print(ind)

    with open(p_log_files[log],'r') as plg:
        pdata = plg.read()
    proc = pdata[pdata.find(' ') : pdata.find('[')]
    pjob = pdata[pdata.find('[')+1 : pdata.find(']')]
#print(pind)

    log_data[ind] = {'lib' : data['library_types'], 'readsn' : data["num_processed"], 'mapn' : data["num_mapped"], 'mapr' : data["percent_mapped"], 'run' : proc, 'job' : pjob }


out_dir = 'logs_' + dirf[dirf.find( 'salmon')+7 : dirf.find('/sra_comp') ] + '.csv'
print(out_dir)
with open(out_dir, 'w') as f:
    f.write("%s,%s,%s,%s,%s,%s,%s\n"%('exp_num','lib_types','reads_processed','reads_mapped','mapping_rate','run','job'))
    for key in log_data.keys():
        f.write("%s,%s,%s,%s,%s,%s,%s\n"%(key,log_data[key]['lib'],log_data[key]['readsn'],log_data[key]['mapn'],log_data[key]['mapr'], log_data[key]['run'], log_data[key]['job']))
