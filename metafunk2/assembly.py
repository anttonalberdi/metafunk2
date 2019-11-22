#!/usr/bin/env python
# -*- coding: utf-8

"""The script to assemble metagenomic reads"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip
import shutil

def assembly(outpath,name,logfilepath,statsfilepath,threads,memory,keep):
    prevdir = "genome_mapping"
    absprevdirr = os.path.join(outpath, name + '.' + prevdir)
    #Create assembly directory if it does not exist
    assembly_dir = "assembly"
    assembly_abs = os.path.join(outpath, name + '.' + assembly_dir)
    #Create genomes directory if it does not exist
    if not os.path.exists(assembly_abs):
        os.makedirs(assembly_abs)
    #Declare input files
    read1in = os.path.join(outpath, name +  '.1.fq')
    read2in = os.path.join(outpath, name +  '.2.fq')

    #Run assembly
    assemblyCmd = 'metaspades.py -1 '+read1in+' -2 '+read2in+' -t '+threads+' -m '+memory+' -k 21,29,39,59,79,99,119 --only-assembler -o '+assembly_abs+''
    subprocess.check_call(assemblyCmd, shell=True)

    #Move reads to working directory
    assembly = os.path.join(assembly_abs, 'contigs.fasta')
    assemblyfinal = os.path.join(outpath, name +  '.fna')
    shutil.copy(assembly, assemblyfinal)

    #Print error to log file if final files are not created
    if not os.path.exists(assemblyfinal):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    There was an error during the assembly. Check error file. \r\n".format(current_time,refgenname))
        logfile.close()

    #If keep is not selected, remove previous directory
    if not keep:
        shutil.rmtree(absprevdirr)
