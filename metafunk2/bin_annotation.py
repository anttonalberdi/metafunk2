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

def bin_annotation(projectpath,threads,memory,logfilepath,threads):
    binannotdir = os.path.join(projectpath, 'merged/binning/annotation')
    if not os.path.exists(binannotdir):
        os.makedirs(binannotdir)
    bindir = os.path.join(projectpath, 'merged/binning')
    databasedir  = '/home/projects/ku-cbd/people/antalb/databases/CAT_prepare_20190719/2019-07-19_CAT_database'
    taxonomydir = '/home/projects/ku-cbd/people/antalb/databases/CAT_prepare_20190719/2019-07-19_taxonomy/'
    proteinfile = os.path.join(projectpath, 'merged/binning/dastool/dastool_proteins.faa')
    outbase = os.path.join(binannotdir, 'annotation')

    batCmd = 'module unload tools gcc/5.1.0 && module load tools anaconda3/4.0.0 diamond/0.9.29 prodigal/2.6.3 cat/5.0.3 && CAT bins -b '+bindir+' -s .fa -d '+databasedir+' -t '+taxonomydir+' -p '+proteinfile+' -o '+outbase+' -n '+threads+' --force '
    subprocess.check_call(batCmd, shell=True)
