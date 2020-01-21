#!/usr/bin/env python
# -*- coding: utf-8

"""The script for contig binning"""
"""Try to upload this into my branch"""
import os
import sys
import random
import argparse
import subprocess
import time
import gzip
import signal

def binning(outpath,name,logfilepath,threads):
    newdir = "binning"
    absnewdir = os.path.join(outpath, name + '.' + newdir)
    if not os.path.exists(absnewdir):
        os.makedirs(absnewdir)

    #Declare input files
    assemblypath = os.path.join(outpath, name + '.fna')
    assemblybampath = os.path.join(outpath,name + '.assembly_mapping', name + '.mapped.bam')

    #########################
    ######## Metabat ########x
    #########################

    metabatdir = os.path.join(absnewdir, 'metabat')
    if not os.path.exists(metabatdir):
        os.makedirs(metabatdir)
    metabatdepthfile = os.path.join(metabatdir, name + '.depth.txt')
    metabatbinbase = os.path.join(metabatdir, name + '.bin')

    #Generate depth file
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Generating metabat depth file from the reads mapped to the assembly \r\n".format(current_time))
    logfile.close()
    metabatdepthfileCmd = 'module unload gcc && module load tools perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth '+metabatdepthfile+' '+assemblybampath+''
    subprocess.check_call(metabatdepthfileCmd, shell=True)

    #Run metabat
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Running metabat binning\r\n".format(current_time))
    logfile.close()
    metabatCmd = 'module unload gcc && module load tools perl/5.20.2 metabat/2.12.1 && metabat2 -i '+assemblypath+' -a '+metabatdepthfile+' -o '+metabatbinbase+' -m 1500 -t '+threads+' --unbinned '
    subprocess.check_call(metabatCmd, shell=True)

    #########################
    ######## Maxbin #########
    #########################

    maxbindir = os.path.join(absnewdir, 'maxbin')
    if not os.path.exists(maxbindir):
        os.makedirs(maxbindir)
    maxbindepthfile = os.path.join(maxbindir, name + '.depth.txt')
    maxbinbase = os.path.join(maxbindir, name + '.bin')

    #Generate depth file
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Generating maxbin depth file from the reads mapped to the assembly \r\n".format(current_time))
    logfile.close()
    maxbindepthfileCmd = 'module unload gcc && module load tools perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth '+maxbindepthfile+' --noIntraDepthVariance '+assemblybampath+''
    subprocess.check_call(maxbindepthfileCmd, shell=True)

    #Run maxbin
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Running maxbin \r\n".format(current_time))
    logfile.close()
    maxbinCmd = 'module unload gcc && module load tools perl/5.20.2 maxbin/2.2.7 fraggenescan/1.31 && run_MaxBin.pl -contig '+assemblypath+' -abund '+maxbindepthfile+' -out '+maxbinbase+' -thread '+threads+''
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
