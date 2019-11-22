#!/usr/bin/env python
# -*- coding: utf-8

"""The script to remove duplicated sequences"""

import os
import sys
import random
import argparse
import subprocess
import gzip

def duplicate_removal(read1,read2,outpath,name,threads,statsfilepath,keep):
    #Create quality_filtering subdirectory
    prevdir = "quality_filtering"
    absprevdirr = os.path.join(outpath, name + '.' + prevdir)
    newdir = "duplicate_removal"
    absnewdir = os.path.join(outpath, name + '.' + newdir)
    if not os.path.exists(absnewdir):
        os.makedirs(absnewdir)

    #Declare input and output files
    read1in = os.path.join(absprevdirr, name +  '.1.fq')
    read2in = os.path.join(absprevdirr, name +  '.2.fq')
    read1tempout = os.path.join(absnewdir, name +  '.1.temp.fq')
    read2tempout = os.path.join(absnewdir, name +  '.2.temp.fq')
    read1out = os.path.join(absnewdir, name +  '.1.fq')
    read2out = os.path.join(absnewdir, name +  '.2.fq')

    #Run seqkit rmdup
    Rmdup1Cmd = 'cat '+read1in+' | seqkit rmdup -s -d bla -o '+read1tempout+''
    subprocess.check_call(Rmdup1Cmd, shell=True)
    Rmdup2Cmd = 'cat '+read2in+' | seqkit rmdup -s -o '+read2tempout+''
    subprocess.check_call(Rmdup2Cmd, shell=True)

    #Repa
    RepCmd = 'repair.sh in='+read1tempout+' in2='+read2tempout+' out='+read1out+' out2='+read2out+' overwrite=t'
    subprocess.check_call(RepCmd, shell=True)

    #Remove temporal files
    os.remove(read1tempout)
    os.remove(read2tempout)

    #Get stats
    reads = 0
    bases = 0
    with open(read1out, 'rb') as read:
        for id in read:
            seq = next(read)
            reads += 1
            bases += len(seq.strip())
            next(read)
            next(read)
    statsfile=open(statsfilepath,"a+")
    statsfile.write("Reads after duplicate removal\t{0}\r\nBases after duplicate removal\t{1}".format(reads,bases))
    statsfile.close()

    #If keep is not selected, remove previous directory
    if not keep:
        os.rmdir(prevdir)
