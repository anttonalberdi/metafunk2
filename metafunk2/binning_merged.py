#!/usr/bin/env python
# -*- coding: utf-8

"""The script for contig binning from reassembly"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip

def binning_merged(projectname,projectpath,threads,memory,logfilepath):
    binningdir = "binning"
    binningdir_abs = os.path.join(projectpath, 'merged',binningdir)
    if not os.path.exists(binningdir_abs):
        os.makedirs(binningdir_abs)

    #Declare input files
    reassemblypath = os.path.join(projectpath, 'merged', 'reassembly.fna')
    reassemblybampaths = os.path.join(projectpath, 'merged','reassembly_mapping','*.sorted.bam')

    #########################
    ######## Metabat ########x
    #########################

    metabatdir = os.path.join(binningdir_abs, 'metabat')
    if not os.path.exists(metabatdir):
        os.makedirs(metabatdir)
    metabatdepthfile = os.path.join(metabatdir, 'depth.txt')
    metabatbinbase = os.path.join(metabatdir, 'bin')

    #Generate depth file
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Generating metabat depth file from the reads mapped to the reassembly \r\n".format(current_time))
    logfile.close()
    metabatdepthfileCmd = 'module unload gcc && module load perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth '+metabatdepthfile+' '+reassemblybampaths+''
    subprocess.check_call(metabatdepthfileCmd, shell=True)

    #Run metabat
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Running metabat binning\r\n".format(current_time))
    logfile.close()
    metabatCmd = 'module unload gcc && module load perl/5.20.2 metabat/2.12.1 && metabat2 -i '+reassemblypath+' -a '+metabatdepthfile+' -o '+metabatbinbase+' -m 1500 -t '+threads+' --unbinned '
    #subprocess.check_call(metabatCmd, shell=True)

    #########################
    ######## Maxbin ######### ERROR! requires FragGeneScan. Waiting for Computerome admin to install it
    #########################

    maxbindir = os.path.join(binningdir_abs, 'maxbin')
    if not os.path.exists(maxbindir):
        os.makedirs(maxbindir)
    maxbindepthfile = os.path.join(maxbindir, 'depth.txt')

    #Generate depth file
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Generating maxbin depth file from the reads mapped to the assembly \r\n".format(current_time))
    logfile.close()
    maxbindepthfileCmd = 'module unload gcc && module load perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth '+maxbindepthfile+' --noIntraDepthVariance '+assemblybampath+''
    subprocess.check_call(maxbindepthfileCmd, shell=True)

    #Run maxbin
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Running maxbin \r\n".format(current_time))
    logfile.close()
    maxbinCmd = 'module load maxbin/2.2.7 perl/5.20.2 fraggenescan/1.31 && run_MaxBin.pl -contig '+assemblypath+' -abund '+maxbindepthfile+' -out '+maxbindir+' -thread '+threads+''
    subprocess.check_call(maxbinCmd, shell=True)

    #######################
    ######## MyCC ######### 2019/11/23 - yelding an error: ValueError: invalid literal for int() with base 10: 'Traceback (most recent call last):\n  File "/services/tools/mycc/20170301/GetThr.py", line 22, in <module>\n    print sorted(dlist,reverse=True)[thr]\nIndexError: list index out of range'
    #######################

    #myccdir = os.path.join(absnewdir, 'mycc')
    #if not os.path.exists(myccdir):
    #    os.makedirs(myccdir)
    #myccdepthfile = os.path.join(myccdir, name + '.depth.txt')

    #Generate depth file
    #logfile=open(logfilepath,"a+")
    #current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    #logfile.write("{0} |    Generating mycc depth file from the reads mapped to the assembly \r\n".format(current_time))
    #logfile.close()
    #metabatdepthfileCmd = 'jgi_summarize_bam_contig_depths --outputDepth '+metabatdepthfile+' '+assemblybampath+''
    #subprocess.check_call(metabatdepthfileCmd, shell=True)

    #Run MyCC
    #logfile=open(logfilepath,"a+")
    #current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    #logfile.write("{0} |    Running MyCC binning\r\n".format(current_time))
    #logfile.close()
    #myccCmd = 'MyCC.py '+assemblypath+' -a '+metabatdepthfile+' '
    #subprocess.check_call(myccCmd, shell=True)

def bin_refinement(projectname,projectpath,threads,memory,logfilepath):

    dastoolDependencies = 'module load gcc/5.4.0 intel/perflibs/2018 R/3.6.1 ruby/2.6.3 pullseq/1.0.2 perl/5.24.0 ncbi-blast/2.6.0+ das_tool/1.1.1'
    dastoolCmd = ''+dastoolDependencies+' && DAS_Tool '
    #https://github.com/cmks/DAS_Tool


    binningdir = "binnin
