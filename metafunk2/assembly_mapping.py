#!/usr/bin/env python
# -*- coding: utf-8

"""The script for contig binning"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip

def assembly_mapping(outpath,name,logfilepath,threads):
    newdir = "assembly_mapping"
    absnewdir = os.path.join(outpath, name + '.' + newdir)
    if not os.path.exists(absnewdir):
        os.makedirs(absnewdir)

    assemblypath = os.path.join(outpath, name + '.fna')
    #Print error and quit if the assembly file does not exist
    if not os.path.exists(assemblypath):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} | ERROR! Assembly file of sample {1} does not exist. \r\n".format(current_time,name))
        logfile.close()
        os.kill(os.getpid(), signal.SIGSTOP)

    assemblybampath = os.path.join(absnewdir, name + '.mapped.bam')
    read1in = os.path.join(outpath, name + '.1.fq')
    read2in = os.path.join(outpath, name + '.2.fq')

    #Index assembly
    assemblyfai = os.path.join(assemblypath + '.fai')
    if not os.path.exists(assemblyfai):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Indexing metagenomic assembly \r\n".format(current_time))
        logfile.close()
        samtoolsindexCmd = 'module load tools samtools/1.9 && samtools faidx '+assemblypath+''
        bwaindexCmd = 'module load tools bwa/0.7.15 && bwa index '+assemblypath+''
        subprocess.check_call(samtoolsindexCmd, shell=True)
        subprocess.check_call(bwaindexCmd, shell=True)
    else:
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Metagenomic assembly is already indexed\r\n".format(current_time))
        logfile.close()

    #Declare mapping commands
    mapCmd = 'module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+threads+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+assemblypath+' '+read1in+' '+read2in+' | samtools view -T '+assemblypath+' -b - | samtools sort -T '+assemblypath+' - > '+assemblybampath+''

    #Mapping to genome
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Mapping reads to assembly \r\n".format(current_time))
    logfile.close()
    subprocess.check_call(mapCmd, shell=True)
