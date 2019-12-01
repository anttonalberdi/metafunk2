#!/usr/bin/env python
# -*- coding: utf-8

"""The script to quality-filter the sequencing reads"""

import os
import sys
import random
import argparse
import subprocess
import gzip

def quality_filtering(read1,read2,outpath,name,threads,statsfilepath,logfilepath):
    #Create quality_filtering subdirectory
    subdir = "quality_filtering"
    absdir = os.path.join(outpath, name + '.' + subdir)
    if not os.path.exists(absdir):
        os.makedirs(absdir)

    #Get initial stats
    reads = 0
    bases = 0
    #If gzipped
    if read1.endswith('.gz'):
        with gzip.open(read1, 'rb') as read:
            for id in read:
                seq = next(read)
                reads += 1
                bases += len(seq.strip())
                next(read)
                next(read)
    else:
        with open(read1, 'rb') as read:
            for id in read:
                seq = next(read)
                reads += 1
                bases += len(seq.strip())
                next(read)
                next(read)
    statsfile=open(statsfilepath,"a+")
    statsfile.write("Input reads\t{0} ({1} bases)\r\n".format(reads,bases))
    statsfile.close()

    #Run Adapterremoval
    ARCmd = 'module unload gcc tools ngs && module load tools gcc/5.4.0 AdapterRemoval/2.1.3 && AdapterRemoval --file1 '+read1+' --file2 '+read2+' --basename '+outpath+'/'+name+'.quality_filtering/'+name+' --minquality 30 --trimqualities --trimns --maxns 5 --threads '+threads+''
    subprocess.check_call(ARCmd, shell=True)

    #Modify output names
    read1in = os.path.join(absdir, name + '.pair1' + '.truncated')
    read2in = os.path.join(absdir, name + '.pair2' + '.truncated')
    read1out = os.path.join(absdir, name + '.1' + '.fq')
    read2out = os.path.join(absdir, name + '.2' + '.fq')
    os.rename(read1in,read1out)
    os.rename(read2in,read2out)

    #Remove unnecesary files
    discarded = os.path.join(absdir, name +  '.discarded')
    os.remove(discarded)
    settings = os.path.join(absdir, name +  '.settings')
    os.remove(settings)
    singleton = os.path.join(absdir, name + '.singleton' + '.truncated')
    os.remove(singleton)

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
    statsfile.write("Quality filtered reads\t{0} ({1} bases)\r\n".format(reads,bases))
    statsfile.close()

    #Print stats to logfile
    logfile=open(logfilepath,"a+")
    logfile.write("                     {0} reads ({1} bases) were kept after quality filtering\r\n".format(reads,bases))
    logfile.close()

    #Doublecheck if everything is ok
    if ( os.stat(read1out).st_size == 0 or  os.stat(read1out).st_size == 0):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} | ERROR! Metafunk2 has stopped due to an error. Check error file \r\n".format(current_time))
        logfile.close()
        os.kill(os.getpid(), signal.SIGSTOP)
