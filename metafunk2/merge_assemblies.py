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
    #Create merged and merged/assembly directories if do not exist
    merged_dir = "merged"
    merged_abs = os.path.join(projectpath, merged_dir)
    if not os.path.exists(merged_abs):
        os.makedirs(merged_abs)
    assembly_dir = "assembly"
    assembly_abs = os.path.join(merged_abs, assembly_dir)
    if not os.path.exists(assembly_abs):
        os.makedirs(assembly_abs)

    assembliespath = os.path.join(projectpath,'*.assembly', 'contigs.fasta')
    mergedassembliespath = os.path.join(assembly_abs, 'assemblies.fna')
    mergedassembliesbase = os.path.join(assembly_abs, 'assemblies')
    nrassembliespath = os.path.join(assembly_abs, 'assemblies.nr.fna')
    afgassembliespath = os.path.join(assembly_abs, 'assemblies.afg')

    #Concatenate assemblies
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Concatenating assemblies \r\n".format(current_time))
    logfile.close()
    concCmd = 'cat '+assembliespath+' > '+mergedassembliespath+''
    subprocess.check_call(concCmd, shell=True)

    #Removing redundant contigs
    #cdhitCmd = 'cd-hit -i '+mergedassembliespath+' -o '+nrassembliespath+' -T '+threads+' -M 0 -c 0.99 -d 100 -aS 0.9'
    #subprocess.check_call(cdhitCmd, shell=True)

    #Modify merged assembly to afg format
    toamosCmd = 'toAmos -s '+mergedassembliespath+' -o '+afgassembliespath+''
    subprocess.check_call(toamosCmd, shell=True)

    #Reassemble assemblies - 20191122 show-coords path issue
    minimusCmd = 'minimus2 '+mergedassembliesbase+' -D OVERLAP=100 -D MINID=95 -D THREADS='+threads+''
    #subprocess.check_call(minimusCmd, shell=True)
