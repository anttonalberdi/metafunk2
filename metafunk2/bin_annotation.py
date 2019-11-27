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
import glob

def bin_annotation(projectpath,threads,memory,logfilepath):
    binannotdir = os.path.join(projectpath, 'merged/binning/annotation')
    if not os.path.exists(binannotdir):
        os.makedirs(binannotdir)
    bindir = os.path.join(projectpath, 'merged/binning')
    databasedir  = '/home/projects/ku-cbd/people/antalb/databases/CAT_prepare_20190719/2019-07-19_CAT_database'
    taxonomydir = '/home/projects/ku-cbd/people/antalb/databases/CAT_prepare_20190719/2019-07-19_taxonomy/'
    outbase = os.path.join(binannotdir, 'annotation')

    ############
    # Taxonomic annotation
    ############

    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |        Running taxonomic annotation using CAT. \r\n".format(current_time))
    logfile.close()

    #Run annotation
    batCmd = 'module unload tools gcc/5.1.0 && module load tools anaconda3/4.0.0 diamond/0.9.29 prodigal/2.6.3 cat/5.0.3 && CAT bins -b '+bindir+' -s .fa -d '+databasedir+' -t '+taxonomydir+' -o '+outbase+' -n '+threads+' --force '
    #subprocess.check_call(batCmd, shell=True)

    #Add names to taxids
    orffile = os.path.join(projectpath, 'merged/binning/annotation/annotation.bin2classification.txt')
    taxfile = os.path.join(projectpath, 'merged/binning/bin.annotation.txt')
    taxonomydir = '/home/projects/ku-cbd/people/antalb/databases/CAT_prepare_20190719/2019-07-19_taxonomy/'
    taxCmd = 'module unload tools gcc/5.1.0 && module load tools anaconda3/4.0.0 diamond/0.9.29 prodigal/2.6.3 cat/5.0.3 && CAT add_names -i '+orffile+' -o '+taxfile+' -t '+taxonomydir+' --only_official'
    #subprocess.check_call(taxCmd, shell=True)

    #Sumarize
    summaryfile = os.path.join(projectpath, 'merged/binning/bin.annotation.summary.txt')
    summaryCmd = 'module unload tools gcc/5.1.0 && module load tools anaconda3/4.0.0 diamond/0.9.29 prodigal/2.6.3 cat/5.0.3 && CAT summarise -i '+taxfile+' -o '+summaryfile+''
    #subprocess.check_call(summaryCmd, shell=True)

    ############
    # Functional annotation
    ############

    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |        Running functional annotation using prokka. \r\n".format(current_time))
    logfile.close()

    binlist = glob.glob(os.path.join(projectpath, 'merged/binning/*.fa'))
    bincount = len(binlist)
    for i in range(bincount):
        bin = binlist[i]
        binname = os.path.splitext(os.path.basename(bin))[0]+''
        outdir = os.path.join(projectpath, 'merged/binning/annotation/' + binname)
        prokkaCmd = 'module load tools signalp/4.1c perl/5.24.0 prokka/1.14.0 && prokka --outdir '+outdir+' --prefix '+binname+' --cpus '+threads+' --force '+bin+''
        subprocess.check_call(prokkaCmd, shell=True)
