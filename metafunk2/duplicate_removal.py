#!/usr/bin/env python
# -*- coding: utf-8

"""The script to remove duplicated sequences"""

import os
import sys
import random
import argparse
import subprocess
import gzip
import shutil

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

    #Print stats to stats file
    statsfile=open(statsfilepath,"a+")
    statsfile.write("Dereplicated reads\t{0} ({1} bases)\r\n".format(reads,bases))
    statsfile.close()

    #Print stats to logfile
    logfile=open(logfilepath,"a+")
    logfile.write("                 {0} reads ({1} bases) were kept after duplicate removal\r\n".format(reads,bases))
    logfile.close()

    #If keep is not selected, remove previous directory
    if not keep:
        if os.path.exists(absprevdirr):
            shutil.rmtree(absprevdirr)

    #Doublecheck everything is ok
    if ( not os.path.exists(read1out) or not os.path.exists(read2out) ):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} | ERROR! Metafunk2 has stopped due to an error. Check error file \r\n".format(current_time))
        logfile.close()
        os.kill(os.getpid(), signal.SIGSTOP)
