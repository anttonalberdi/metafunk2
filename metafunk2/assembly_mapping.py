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
import signal

def assembly_mapping(outpath,name,logfilepath,threads):
    newdir = "assembly_mapping"
    absnewdir = os.path.join(outpath, name + '.' + newdir)
    if not os.path.exists(absnewdir):
        os.makedirs(absnewdir)

    #Print to log
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2 has started mapping the reads back to the metagenomic assembly \r\n".format(current_time))
    logfile.close()

    ###########
    # 1) DECLARE INPUT AND OUTPUT FILES AND CHECK IF THEY EXIST
    ###########
    # The assembly requires two input elements and one output element to be declared
    # 1) Assembly file, 2) Forward and reverse read files to be mapped to the assembly,
    # and 3) the output bam file.

    #Declare the assembly file
    assemblypath = os.path.join(outpath, name + '.fna')

    #Print error and kill process if the assembly file does not exist
    if not os.path.exists(assemblypath):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} | ERROR! Assembly file of sample {1} does not exist. \r\n".format(current_time,name))
        logfile.close()
        os.kill(os.getpid(), signal.SIGSTOP)

    #Declare the read files
    read1in = os.path.join(outpath, name + '.1.fq')
    read2in = os.path.join(outpath, name + '.2.fq')
    #Print error and kill process if the read files do not exist
    if (not os.path.exists(read1in) or not os.path.exists(read2in)):
        ## Print to log
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} | ERROR! Read files of sample {1} do not exist. \r\n".format(current_time,name))
        logfile.close()
        ## Kill the job
        os.kill(os.getpid(), signal.SIGSTOP)
        ##

    #Declare the output bam file
    assemblybampath = os.path.join(absnewdir, name + '.mapped.bam')

    ###########
    # 2) CHECK WHETHER THE ASSEMBLY IS INDEXED AND OTHERWISE INDEX IT
    ###########
    # The assembly needs to be indexed in order the reads to be mapped to it.

    #Run the indexing only if the .fai file does not exist
    assemblyfai = os.path.join(assemblypath + '.fai')
    if not os.path.exists(assemblyfai):
        ## Print to log
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Indexing metagenomic assembly \r\n".format(current_time))
        logfile.close()
        ##
        samtoolsindexCmd = 'module load tools samtools/1.9 && samtools faidx '+assemblypath+''
        bwaindexCmd = 'module load tools bwa/0.7.15 && bwa index '+assemblypath+''
        subprocess.check_call(samtoolsindexCmd, shell=True)
        subprocess.check_call(bwaindexCmd, shell=True)
    else:
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Metagenomic assembly is already indexed\r\n".format(current_time))
        logfile.close()

    ###########
    # 3) MAP THE READS TO THE ASSEMBLY
    ###########

    ## Print to log
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Mapping reads to {1} assembly \r\n".format(current_time,name))
    logfile.close()
    ##

    mapCmd = 'module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+threads+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+assemblypath+' '+read1in+' '+read2in+' | samtools view -T '+assemblypath+' -b - | samtools sort -T '+assemblypath+' - > '+assemblybampath+''
    subprocess.check_call(mapCmd, shell=True)
