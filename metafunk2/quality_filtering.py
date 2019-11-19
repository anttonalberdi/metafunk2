#!/usr/bin/env python
# -*- coding: utf-8

"""The script to quality-filter the sequencing reads"""

import os
import sys
import random
import argparse
import subprocess

def quality_filtering(read1,read2,outpath,name,threads):
    #Create quality_filtering subdirectory
    subdir = "quality_filtering"
    absdir = os.path.join(outpath, name + '.' + subdir)
    if not os.path.exists(absdir):
        os.makedirs(absdir)

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
