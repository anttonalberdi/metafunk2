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
        #Declare genome name and path
        refgenname = refgenlist[i][0]
        refgenpath = os.path.join(outpath,'genomes', refgenname + '.fna')
        refgenfai = os.path.join(outpath,'genomes', refgenname + '.fna.fai')
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

def genome_mapping(refgenlist,outpath,name,logfilepath,threads):
    #Create quality_filtering subdirectory
    prevdir = "duplicate_removal"
    absprevdirr = os.path.join(outpath, name + '.' + prevdir)

    #Declare input files
    read1in = os.path.join(absprevdirr, name +  '.1.fq')
    read2in = os.path.join(absprevdirr, name +  '.2.fq')

    #Iterate across reference genomes
    refgencount = len(refgenlist)
    for i in range(refgencount):
        if i>0:
            read1in = read1out
            read2in = read2out
        #Declare genome name and path
        refgenname = refgenlist[i][0]
        refgenpath = os.path.join(outpath,'genomes', refgenname + '.fna')
        bampath_all = os.path.join(outpath,'genome_mapping', name + '.' + refgenname + '.bam')
        bampath_host = os.path.join(outpath,'genome_mapping', name + '.mappedto.' + refgenname + '.bam')
        bampath_mg = os.path.join(outpath,'genome_mapping', name + '.mg.bam')
        singletonpath = os.path.join(outpath,'genome_mapping', name + '.singleton.fq')
        read1out = os.path.join(outpath,'genome_mapping', name +  '.1.fq')
        read2out = os.path.join(outpath,'genome_mapping', name +  '.2.fq')

        #Declare mapping commands
        mapCmd = 'bwa mem -t '+threads+' '+refgenpath+' '+read1in+' '+read2in+' | samtools view -T '+refgenpath+' -b - > '+bampath_all+''
        hostmapCmd = 'samtools view -T '+refgenpath+' -b -F12 '+bampath_all+' > '+bampath_host+''
        mgmapCmd = 'samtools view -T '+refgenpath+' -b -f12 '+bampath_all+' > '+bampath_mg+''
        mgfqCmd = 'samtools fastq -s '+singletonpath+' -1 '+read1out+' -2 '+read2out+' '+bampath_mg+''

        #Mapping to genome
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Mapping reads to {1} genome \r\n".format(current_time,refgenname))
        logfile.close()
        subprocess.check_call(mapCmd, shell=True)

        #Extracting mapped (genomic) reads
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Extracting reads mapped to {1} genome (genomic reads) \r\n".format(current_time,refgenname))
        logfile.close()
        subprocess.check_call(hostmapCmd, shell=True)

        #Extracting unmapped (metagenomic) reads
        #Extracting mapped (genomic) reads
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Extracting reads not mapped to {1} genome (metagenomic reads) \r\n".format(current_time,refgenname))
        logfile.close()
        #Extract unmapped
        subprocess.check_call(mgmapCmd, shell=True)
        #Convert unmapped to fastq
        subprocess.check_call(mgfqCmd, shell=True)
