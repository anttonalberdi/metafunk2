#!/usr/bin/env python
# -*- coding: utf-8

"""The script to map read to reference genomes"""

import os
import sys
import random
import argparse
import subprocess
from shutil import copyfile

#Copy reference genome to working directory
def copy_genome(refgenomepath,outpath):
    #Create genome_mapping directory if it does not exist
    newdir = "genome_mapping"
    absnewdir = os.path.join(outpath, name + '.' + newdir)
    #Create genome_mapping/genomes directory if it does not exist
    if not os.path.exists(absnewdir):
        os.makedirs(absnewdir)
    genomedir = "genomes"
    absgenomedir = os.path.join(absnewdir, genomedir)
    if not os.path.exists(absgenomedir):
        os.makedirs(absgenomedir)
    #Copy genome file to genome_mapping/genomes directory
    copyfile(refgenomepath, absgenomedir)
    refgen = os.path.basename(refgenomepath)

#Manipulate reference genome
def check_genome(refgenomepath):
    if refgenpath.endswith('.gz'):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("     {0} | Decompressing reference genome {1} \r\n".format(current_time,refgenomepath))
        logfile.close()
        DecompCmd = 'pigz '+refgenpath+''
        subprocess.check_call(DecompCmd, shell=True)
        refgenomepathnew = refgenomepath.replace(".gz", "")
    if refgenpath.endswith('.fasta'):
        original = os.path.join(absdir, name + '.pair1' + '.truncated')
        os.rename(read1in,read1out)
    if refgenpath.endswith('.fa'):
        os.rename(read1in,read1out)

#def index_genome(refgenomepath):

#def genome_mapping(read1,read2,refgenomepath):
