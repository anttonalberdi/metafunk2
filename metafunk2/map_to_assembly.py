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

def map_to_assembly(outpath,name,logfilepath,threads):
    newdir = "binning"
    absnewdir = os.path.join(outpath, name + '.' + newdir)
    if not os.path.exists(absnewdir):
        os.makedirs(absnewdir)

    assemblypath = os.path.join(outpath, name + '.assembly','contigs.fasta')
    assemblybampath = os.path.join(outpath,name + '.genome_mapping', name + '.mapped.bam')
    read1in = os.path.join(outpath, name + '.genome_mapping',name + '.1.fna')
    read2in = os.path.join(outpath, name + '.genome_mapping',name + '.2.fna')
    #Index assembly
    assemblyfai = os.path.join(assemblypath + '.fai')
    if not os.path.exists(assemblyfai):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Indexing metagenomic assembly \r\n".format(current_time))
        logfile.close()
        samtoolsindexCmd = 'samtools faidx '+assemblypath+''
        bwaindexCmd = 'bwa index '+assemblypath+''
        subprocess.check_call(samtoolsindexCmd, shell=True)
        subprocess.check_call(bwaindexCmd, shell=True)
    else:
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Metagenomic assembly is already indexed\r\n".format(current_time))
        logfile.close()

    #Declare mapping commands
    mapCmd = 'bwa mem -t '+threads+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+assemblypath+' '+read1in+' '+read2in+' | samtools sort -T '+assemblypath+' -b - > '+assemblybampath+''

    #Mapping to genome
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Mapping reads to assemblt \r\n".format(current_time,refgenname))
    logfile.close()
    subprocess.check_call(mapCmd, shell=True)
