#!/usr/bin/env python
# -*- coding: utf-8

"""The script for contig binning"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip

def bin_qc(projectpath,threads,memory,logfilepath):
    bindir = os.path.join(projectpath, 'merged/binning/metabat')
    bindirtemp = os.path.join(projectpath, 'merged/binning/metabat.tmp')
    if not os.path.exists(bindirtemp):
        os.makedirs(bindirtemp)
    outputdir = os.path.join(projectpath, 'merged/binning/metabat.checkm')
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    checkmCmd = 'module load anaconda2/4.0.0 && checkm lineage_wf -t '+threads+' --pplacer_threads '+threads+' --tmpdir '+bindirtemp+' '+bindir+' '+outputdir+''
    subprocess.check_call(checkmCmd, shell=True)
