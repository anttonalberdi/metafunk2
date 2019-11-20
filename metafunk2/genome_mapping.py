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
def copy_genome(refgenomepath,outpath,name,logfilepath):
    #Create genome_mapping directory if it does not exist
    newdir = "genome_mapping"
    absnewdir = os.path.join(outpath, name + '.' + newdir)
    #Create genomes directory if it does not exist
    if not os.path.exists(absnewdir):
        os.makedirs(absnewdir)
    genomedir = "genomes"
    absgenomedir = os.path.join(outpath, genomedir)
    if not os.path.exists(absgenomedir):
        os.makedirs(absgenomedir)

    #Copy genome file to genomes directory if it is not already there
    refgen = os.path.basename(refgenomepath)
    newgenomepath = os.path.join(absgenomedir, refgen)
    newgenomepathnoext = newgenomepath.rsplit( ".", 1 )[ 0 ]
    newgenomepathnoext2 = newgenomepath.rsplit( ".", 1 )[ 0 ]
    newgenomepath2 = os.path.join(newgenomepathnoext + '.fna')
    newgenomepath3 = os.path.join(newgenomepathnoext2 + '.fna')
    print(newgenomepath)
    print(newgenomepath2)
    if ( not os.path.exists(newgenomepath) and not os.path.exists(newgenomepath2) and not os.path.exists(newgenomepath3) ):
        #Add to log
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Transferring genome {1} to working directory \r\n".format(current_time,refgenomepath))
        logfile.close()
        copyfile(refgenomepath, newgenomepath)
    else:
        #Add to log
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Genome {1} already exists in working directory \r\n".format(current_time,newgenomepath))
        logfile.close()

    #Replace original genome path with new path
    refgenomepath = newgenomepath

    #Manipulate reference genome
    if refgenomepath.endswith('.gz'):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Decompressing reference genome {1} \r\n".format(current_time,refgenomepath))
        logfile.close()
        DecompCmd = 'pigz -d '+refgenomepath+''
        subprocess.check_call(DecompCmd, shell=True)
        refgenomepath = refgenomepath.replace(".gz", "")
    if refgenomepath.endswith('.fasta'):
        refgenomepathnoext = refgenomepath.rsplit( ".", 1 )[ 0 ]
        newrefgenpath = os.path.join(refgenomepathnoext + '.fna')
        print(refgenomepath)
        print(newrefgenpath)
        os.rename(refgenomepath,newrefgenpath)
    if refgenomepath.endswith('.fa'):
        refgenomepathnoext = refgenomepath.rsplit( ".", 1 )[ 0 ]
        newrefgenpath = os.path.join(refgenomepathnoext + '.fna')
        print(refgenomepath)
        print(newrefgenpath)
        os.rename(refgenomepath,newrefgenpath)

#def index_genome(refgenomepath):

#def genome_mapping(read1,read2,refgenomepath):
