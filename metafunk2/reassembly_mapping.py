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
import glob
import ntpath
from check_software import is_tool

def reassembly_indexing(projectname,projectpath,threads,memory,logfilepath):
    #Create reassembly_mapping directory if it does not exist
    reassembly_mapping_dir = "reassembly_mapping"
    reassembly_mapping_dir_abs = os.path.join(projectpath, 'merged', reassembly_mapping_dir)
    #Create genomes directory if it does not exist
    if not os.path.exists(reassembly_mapping_dir_abs):
        os.makedirs(reassembly_mapping_dir_abs)

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
        bwaindexCmd = 'module load bwa/0.7.15 && bwa index '+reassemblypath+''
        subprocess.check_call(samtoolsindexCmd, shell=True)
        subprocess.check_call(bwaindexCmd, shell=True)
    else:
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Metagenomic reassembly (merged assemblies) is already indexed\r\n".format(current_time))
        logfile.close()

def reassembly_mapping(projectname,projectpath,threads,memory,logfilepath):
    #Map reads from each sample to the assembly
    #Detect samples
    read1 = glob.glob(os.path.join(projectpath, '*.1.fq'))
    read2 = glob.glob(os.path.join(projectpath, '*.2.fq'))

    #Error if 1 and 2 reads of all samples are not present
    if len(read1) != len(read2):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    ERROR! The automatically detected number of forward and reverse reads is not the same.\r\n".format(current_time))
        logfile.close()
        os.kill(os.getpid(), signal.SIGSTOP)

    #Add to log
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Mapping reads of all samples to the reassembly\r\n".format(current_time))
    logfile.close()

    #Declare mapping command
    reassemblypath = os.path.join(projectpath, 'merged', 'reassembly.fna')

    samplecount = len(read1)
    for i in range(samplecount):
        read1in = read1[i]
        read2in = read1[i]
        name = os.path.splitext(os.path.basename(reads1in))[0]+''
        bampath_mapped = os.path.join(reassembly_mapping_dir_abs, name + '.bam')

        #Run mapping
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |            Mapping {1} reads to the reassembly\r\n".format(current_time,name))
        logfile.close()
        mapCmd = 'module load samtools/1.9 bwa/0.7.15 && bwa mem -t '+threads+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+reassemblypath+' '+read1in+' '+read2in+' | samtools view -T '+reassemblypath+' -b -f12 - > '+bampath_mapped+''
        subprocess.check_call(mapCmd, shell=True)
