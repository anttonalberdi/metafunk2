#!/usr/bin/env python
# -*- coding: utf-8

"""The script for  bin annotation"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip

def bin_annotation(projectpath,threads,memory,logfilepath):
    binannotdir = os.path.join(projectpath, 'merged/binning/annotation')
    if not os.path.exists(binannotdir):
        os.makedirs(binannotdir)
    bindir = os.path.join(projectpath, 'merged/binning')
    databasedir  = '/home/projects/ku-cbd/people/antalb/databases/CAT_prepare_20190719/2019-07-19_CAT_database'
    taxonomydir = '/home/projects/ku-cbd/people/antalb/databases/CAT_prepare_20190719/2019-07-19_taxonomy/'
    outbase = os.path.join(binannotdir, 'annotation')

    #Run annotation
    batCmd = 'module unload tools gcc/5.1.0 && module load tools anaconda3/4.0.0 diamond/0.9.29 prodigal/2.6.3 cat/5.0.3 && CAT bins -b '+bindir+' -s .fa -d '+databasedir+' -t '+taxonomydir+' -o '+outbase+' -n '+threads+' --force '
    subprocess.check_call(batCmd, shell=True)

    #Add names to taxids
    orffile = os.path.join(projectpath, 'merged/binning/annotation/annotation.bin2classification.txt')
    taxfile = os.path.join(projectpath, 'merged/binning/bin.annotation.txt')
    taxonomydir = '/home/projects/ku-cbd/people/antalb/databases/CAT_prepare_20190719/2019-07-19_taxonomy/'
    taxCmd = 'module unload tools gcc/5.1.0 && module load tools anaconda3/4.0.0 diamond/0.9.29 prodigal/2.6.3 cat/5.0.3 && CAT add_names -i '+orffile+' -o '+taxfile+' -t '+taxonomydir+' --only_official'
    subprocess.check_call(taxCmd, shell=True)

    #Sumarize
    summaryfile = os.path.join(projectpath, 'merged/binning/bin.annotation.summary.txt')
    taxCmd = 'module unload tools gcc/5.1.0 && module load tools anaconda3/4.0.0 diamond/0.9.29 prodigal/2.6.3 cat/5.0.3 && CAT CAT summarise -i '+taxfile+' -o '+summaryfile+''
