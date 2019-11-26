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

#WORKING ON IT
def bin_annotation(projectpath,threads,memory,logfilepath):
    bindir = os.path.join(projectpath, 'merged/binning/')
    bindirtemp = os.path.join(projectpath, 'merged/binning/tmp')
    if not os.path.exists(bindirtemp):
        os.makedirs(bindirtemp)
    outputdir = os.path.join(projectpath, 'merged/binning/checkm')
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    #module unload tools gcc/5.1.0 && module load tools anaconda3/4.0.0 diamond/0.9.29 prodigal/2.6.3 cat/5.0.3 && CAT bins -b

    checkmCmd = 'module load tools anaconda2/4.0.0 hmmer/3.2.1 prodigal/2.6.3 pplacer/1.1.alpha19 && checkm lineage_wf -t '+threads+' --pplacer_threads '+threads+' --tmpdir '+bindirtemp+' -x fa '+bindir+' '+outputdir+''
    subprocess.check_call(checkmCmd, shell=True)
