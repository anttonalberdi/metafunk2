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
from shutil import movefile

def assembly(outpath,name,logfilepath,statsfilepath,threads,memory):
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

    #Move reads to parent folder
    assembly = os.path.join(assembly_dir, 'contigs.fasta')
    assemblyfinal = os.path.join(outpath, name +  '.fna')
    shutil.move(read1out, assemblyfinal)
