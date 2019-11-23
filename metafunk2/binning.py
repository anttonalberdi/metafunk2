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

    #Generate depth file
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Generating metabat depth file from the reads mapped to the assembly \r\n".format(current_time))
    logfile.close()
    metabatdepthfileCmd = 'jgi_summarize_bam_contig_depths --outputDepth '+metabatdepthfile+' '+assemblybampath+''
    subprocess.check_call(metabatdepthfileCmd, shell=True)

    #Run metabat
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Running metabat \r\n".format(current_time))
    logfile.close()
    metabatCmd = 'metabat2 -i '+assemblypath+' -a '+metabatdepthfile+' -o '+metabatdir+' -m 1500 -t '+threads+' --unbinned '
    subprocess.check_call(metabatCmd, shell=True)

    #########################
    ######## Maxbin #########
    #########################

    #maxbindir = os.path.join(absnewdir, 'metabat')
    #if not os.path.exists(maxbindir):
#        os.makedirs(maxbindir)
    #maxbindepthfile = os.path.join(maxbindir, name + '.depth.txt')

    #Generate depth file
    #logfile=open(logfilepath,"a+")
    #current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    #logfile.write("{0} |    Generating maxbin depth file from the reads mapped to the assembly \r\n".format(current_time))
    #logfile.close()
    #metabatdepthfileCmd = 'jgi_summarize_bam_contig_depths --outputDepth '+maxbindepthfile+' --noIntraDepthVariance '+assemblybampath+''
    #subprocess.check_call(metabatdepthfileCmd, shell=True)

    #Run maxbin
    #logfile=open(logfilepath,"a+")
    #current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    #logfile.write("{0} |    Running maxbin \r\n".format(current_time))
    #logfile.close()
    #metabatCmd = 'metabat2 -i '+assemblypath+' -a '+metabatdepthfile+' -o '+metabatdir+' -m 1500 -t '+threads+' --unbinned '
    #subprocess.check_call(metabatCmd, shell=True)
