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
    bindir = os.path.join(projectpath, 'merged/binning/')
    bindirtemp = os.path.join(projectpath, 'merged/binning/tmp')
    if not os.path.exists(bindirtemp):
        os.makedirs(bindirtemp)
    outputdir = os.path.join(projectpath, 'merged/binning/checkm')
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    checkmCmd = 'module unload gcc/5.1.0 && module load tools anaconda2/4.0.0 hmmer/3.2.1 prodigal/2.6.3 pplacer/1.1.alpha19 && checkm lineage_wf -t '+threads+' --pplacer_threads '+threads+' --tmpdir '+bindirtemp+' -x fa '+bindir+' '+outputdir+''
    subprocess.check_call(checkmCmd, shell=True)


    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |            Generating plots of completeness, contamination, and strain heterogeneity. \r\n".format(current_time))
    logfile.close()
    plotsdir = os.path.join(outputdir, 'plots')
    plotsCmd = 'module unload gcc/5.1.0 && module load tools anaconda2/4.0.0 hmmer/3.2.1 prodigal/2.6.3 pplacer/1.1.alpha19 && checkm bin_qa_plot -x fa '+plotsdir+' '+bindir+' '+outputdir+''
    subprocess.check_call(plotsCmd, shell=True)
