#!/usr/bin/env python
# -*- coding: utf-8

"""The script to quality-filter the sequencing reads"""

import os
import sys
import random
import argparse
import subprocess

def quality_filtering(read1,read2,outpath,name,threads,statsfilepath):
    #Create quality_filtering subdirectory
    subdir = "quality_filtering"
    absdir = os.path.join(outpath, name + '.' + subdir)
    if not os.path.exists(absdir):
        os.makedirs(absdir)

    #Get initial stats
    reads = 0
    bases = 0
    with gzip.open(read1, 'rb') as read:
        for id in read:
            seq = next(read)
            reads += 1
            bases += len(seq.strip())
            next(read)
            next(read)
    statsfile=open(statsfilepath,"a+")
    statsfile.write("Initial reads\t{0}\r\nInitial bases\t{1}".format(reads,bases))
    statsfile.close()

    #Run Adapterremoval
    ARCmd = 'AdapterRemoval --file1 '+read1+' --file2 '+read2+' --basename '+outpath+'/'+name+'.quality_filtering/'+name+' --minquality 30 --trimqualities --trimns --maxns 5 --threads '+threads+''
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
    with gzip.open(read1out, 'rb') as read:
        for id in read:
            seq = next(read)
            reads += 1
            bases += len(seq.strip())
            next(read)
            next(read)
    statsfile=open(statsfilepath,"a+")
    statsfile.write("Reads after quality filtering\t{0}\r\nBases after quality filtering\t{1}".format(reads,bases))
    statsfile.close()
