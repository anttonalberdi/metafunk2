#!/usr/bin/env python
# -*- coding: utf-8

"""The script to map read to reference genomes"""

import os
import sys
import random
import argparse
import subprocess
import time
from shutil import copyfile

#Copy reference genome to working directory
def copy_genome(refgenlist,outpath,name,logfilepath):
    #Create genome_mapping directory if it does not exist
    genome_mapping_dir = "genome_mapping"
    genome_mapping_dir_abs = os.path.join(outpath, name + '.' + genome_mapping_dir)
    #Create genomes directory if it does not exist
    if not os.path.exists(genome_mapping_dir_abs):
        os.makedirs(genome_mapping_dir_abs)
    genomes_dir = "genomes"
    genomes_dir_abs = os.path.join(outpath, genomes_dir)
    if not os.path.exists(genomes_dir_abs):
        os.makedirs(genomes_dir_abs)

    #Copy genome file(s) to genomes directory if it is not already there
    refgencount = len(refgenlist)
    for i in range(refgencount):
        #Declare genome name
        refgenname = refgenlist[i][0]
        #Declare original genome path
        refgenoriginalpath = refgenlist[i][1]
        #Declare new genome path and copy file
        if refgenoriginalpath.endswith('.gz'):
            refgenpath = os.path.join(genomes_dir_abs, refgenname + '.fna.gz')
            refgenpathuncomp = refgenpath.replace(".gz", "")
            if ( not os.path.exists(refgenpath) and not os.path.exists(refgenpathuncomp) ):
                logfile=open(logfilepath,"a+")
                current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
                logfile.write("{0} |    Transferring {1} genome to working directory \r\n".format(current_time,refgenname))
                logfile.close()
                copyfile(refgenoriginalpath, refgenpath)
            else:
                logfile=open(logfilepath,"a+")
                current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
                logfile.write("{0} |    {1} genome is already in the working directory \r\n".format(current_time,refgenname))
                logfile.close()
        if ( refgenoriginalpath.endswith('.fasta') or refgenoriginalpath.endswith('.fa') or refgenoriginalpath.endswith('.fna') ):
            refgenpath = os.path.join(genomes_dir_abs, refgenname + '.fna')
            if not os.path.exists(refgenpath):
                logfile=open(logfilepath,"a+")
                current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
                logfile.write("{0} |    Transferring {1} genome to working directory \r\n".format(current_time,refgenname))
                logfile.close()
                copyfile(refgenoriginalpath, refgenpath)
            else:
                logfile=open(logfilepath,"a+")
                current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
                logfile.write("{0} |    {1} genome is already in the working directory \r\n".format(current_time,refgenname))
                logfile.close()
        print(refgenname)
        print(refgenoriginalpath)
        print(refgenpath)

        #Manipulate reference genome
        if refgenpath.endswith('.gz'):
            refgenpathuncomp = refgenpath.replace(".gz", "")
            if not os.path.exists(refgenpathuncomp):
                logfile=open(logfilepath,"a+")
                current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
                logfile.write("{0} |    Decompressing {1} genome \r\n".format(current_time,refgenname))
                logfile.close()
                DecompCmd = 'pigz -d '+refgenpath+''
                subprocess.check_call(DecompCmd, shell=True)

#Index reference genomes if not already indexed
def index_genome(refgenlist,outpath,name,logfilepath):
    refgencount = len(refgenlist)
    for i in range(refgencount):
        #Declare genome name
        refgenname = refgenlist[i][0]
        refgenpath = os.path.join(outpath,'genomes', refgenname + '.fna')
        refgenfai = os.path.join(outpath,'genomes', refgenname + '.fai')
        if not os.path.exists(refgenfai):
            logfile=open(logfilepath,"a+")
            current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
            logfile.write("{0} |    Indexing {1} genome \r\n".format(current_time,refgenname))
            logfile.close()
            samtoolsindexCmd = 'samtools faidx '+refgenpath+''
            bwaindexCmd = 'bwa index '+refgenpath+''
            subprocess.check_call(samtoolsindexCmd, shell=True)
            subprocess.check_call(bwaindexCmd, shell=True)
        else:
            logfile=open(logfilepath,"a+")
            current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
            logfile.write("{0} |    {1} genome is already indexed \r\n".format(current_time,refgenname))
            logfile.close()



#def genome_mapping(read1,read2,refgenomepath):
