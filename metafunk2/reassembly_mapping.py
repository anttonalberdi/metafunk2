#!/usr/bin/env python
# -*- coding: utf-8

"""The script for mapping reads to the merged assembly"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip
from check_software import is_tool

def reassembly_indexing(projectname,projectpath,threads,memory,logfilepath):

    newdir = "reassembly_mapping"
    reassemblypath = os.path.join(projectpath, 'merged', 'reassembly.fna')
    if not os.path.exists(reassemblypath):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Merged assembly (reassembly) file does not exist\r\n".format(current_time))
        logfile.close()
        #Kill process if file does not exist
        os.kill(os.getpid(), signal.SIGSTOP)

    #Index reassembly
    reassemblyfai = os.path.join(reassemblypath + '.fai')
    if not os.path.exists(reassemblyfai):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Indexing metagenomic reassembly (merged assemblies) \r\n".format(current_time))
        logfile.close()
        samtoolsindexCmd = 'module load samtools/1.9 && samtools faidx '+reassemblypath+''
        bwaindexCmd = 'bwa/0.7.15 && bwa index '+reassemblypath+''
        subprocess.check_call(samtoolsindexCmd, shell=True)
        subprocess.check_call(bwaindexCmd, shell=True)
    else:
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Metagenomic reassembly (merged assemblies) is already indexed\r\n".format(current_time))
        logfile.close()

def reassembly_mapping(projectname,projectpath,threads,memory,logfilepath):
    newdir = "reassembly_mapping"
    absnewdir = os.path.join(projectpath, 'merged', newdir)
    if not os.path.exists(absnewdir):
        os.makedirs(absnewdir)

    assembliespath = os.path.join(projectpath, 'merged','assemblies.fna')
    assembliesbampath = os.path.join(absnewdir, name + '.mapped.bam')

    #Index merged assembly
    assembliesfai = os.path.join(assembliespath + '.fai')
    if not os.path.exists(assembliesfai):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Indexing merged assembly \r\n".format(current_time))
        logfile.close()
        samtoolsindexCmd = 'samtools faidx '+assembliespath+''
        bwaindexCmd = 'bwa index '+assembliespath+''
        subprocess.check_call(samtoolsindexCmd, shell=True)
        subprocess.check_call(bwaindexCmd, shell=True)
    else:
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Merged assembly is already indexed\r\n".format(current_time))
        logfile.close()

    #Map reads from each sample to the assembly
    #Detect samples

    #Declare mapping commands
    mapCmd = 'bwa mem -t '+threads+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+assemblypath+' '+read1in+' '+read2in+' | samtools view -T '+assemblypath+' -b - | samtools sort -T '+assemblypath+' - > '+assemblybampath+''

    #Mapping to genome
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Mapping reads to assembly \r\n".format(current_time))
    logfile.close()
    subprocess.check_call(mapCmd, shell=True)
