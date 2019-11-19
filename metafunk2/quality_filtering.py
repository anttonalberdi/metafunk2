#!/usr/bin/env python
# -*- coding: utf-8

"""The script to quality-filter the sequencing reads"""

import os
import sys
import random
import argparse
import subprocess

def quality_filtering(read1,read2,outpath,name,threads):
    subdir = "quality_filtering"
    absdir = os.path.join(outpath, subdir)
    if not os.path.exists(absdir):
        os.makedirs(absdir)
    myCmd = 'AdapterRemoval --file1 '+read1+' --file2 '+read2+' --basename '+outpath+'/quality_filtering/'+name+' --minquality 30 --trimqualities --trimns --maxns 5 --threads '+threads+''
    FNULL = open(os.devnull, 'w')
    subprocess.Popen(myCmd, shell=True)
