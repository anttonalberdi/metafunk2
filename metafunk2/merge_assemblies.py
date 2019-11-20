#!/usr/bin/env python
# -*- coding: utf-8

"""The script for merging assemblies"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip

def merge_assemblies(projectname,projectpath,threads,memory,logfilepath):
    #Create merged directory if it does not exist
    merged_dir = "merged"
    merged_abs = os.path.join(projectpath, merged_dir)
    if not os.path.exists(merged_abs):
        os.makedirs(merged_abs)

    assembliespath = os.path.join(projectpath,'*.assembly', 'contigs.fasta')
    mergedassembliespath = os.path.join(merged_abs, 'assemblies.fna')
    nrassembliespath = os.path.join(merged_abs, 'assemblies.nr.fna')

    #Concatenate assemblies
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Concatenating assemblies \r\n".format(current_time))
    logfile.close()
    concCmd = 'cat '+assembliespath+' > '+mergedassembliespath+''
    subprocess.check_call(concCmd, shell=True)

    #Removing redundant contigs
    cdhitCmd = 'cd-hit -i '+mergedassembliespath+' -o '+nrassembliespath+' -T '+threads+' -M 0 -c 0.99 -d 100 -aS 0.9'
    subprocess.check_call(cdhitCmd, shell=True)
