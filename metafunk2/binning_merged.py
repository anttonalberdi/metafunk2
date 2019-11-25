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
import glob

def binning_merged(projectname,projectpath,threads,memory,logfilepath):
    binningdir = "binning"
    binningdir_abs = os.path.join(projectpath, 'merged',binningdir)
    if not os.path.exists(binningdir_abs):
        os.makedirs(binningdir_abs)

    #Declare input files
    reassemblypath = os.path.join(projectpath, 'merged', 'reassembly.fna')
    reassemblybampaths = os.path.join(projectpath, 'merged','reassembly_mapping','*.sorted.bam')

    #########################
    ######## Metabat ########
    #########################

    metabatdir = os.path.join(binningdir_abs, 'metabat')
    if not os.path.exists(metabatdir):
        os.makedirs(metabatdir)
    metabatdepthfile = os.path.join(metabatdir, 'depth.txt')
    metabatbinbase = os.path.join(metabatdir, 'metabat')

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
    metabatCmd = 'module unload gcc && module load perl/5.20.2 metabat/2.12.1 && metabat2 -i '+reassemblypath+' -a '+metabatdepthfile+' -o '+metabatbinbase+' -m 1500 -t '+threads+''
    subprocess.check_call(metabatCmd, shell=True)

    #Create contig to bin table
    bintablefile = os.path.join(binningdir_abs, 'bins_metabat.txt')
    bintable=open(bintablefile,"a+")
    metabatdir = os.path.join(metabatbinbase + '.*.fa')
    binlist = glob.glob(metabatdir)
    for bin in binlist:
        binname = os.path.splitext(os.path.basename(bin))[0]+''
        with open(bin, 'r') as binfile:
           for line in binfile:
                if line.startswith('>'):
                    contig = line.strip()
                    contig = contig.replace(">", "")
                    bintable.write("{0}\t{1}\r\n".format(contig,binname))
    bintable.close()

    #########################
    ######## Maxbin #########
    #########################

    maxbindir = os.path.join(binningdir_abs, 'maxbin')
    if not os.path.exists(maxbindir):
        os.makedirs(maxbindir)
    maxbindepthfile = os.path.join(maxbindir, 'depth.txt')
    maxbinbase = os.path.join(maxbindir, 'maxbin')

    #Generate depth file
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Generating maxbin depth file from the reads mapped to the assembly \r\n".format(current_time))
    logfile.close()
    maxbindepthfileCmd = 'module unload gcc && module load perl/5.20.2 metabat/2.12.1 && jgi_summarize_bam_contig_depths --outputDepth '+maxbindepthfile+' --noIntraDepthVariance '+reassemblybampaths+''
    subprocess.check_call(maxbindepthfileCmd, shell=True)

    #Run maxbin
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Running maxbin \r\n".format(current_time))
    logfile.close()
    maxbinCmd = 'module load perl/5.20.2 maxbin/2.2.7 fraggenescan/1.31 && run_MaxBin.pl -contig '+reassemblypath+' -abund '+maxbindepthfile+' -out '+maxbinbase+' -thread '+threads+''
    subprocess.check_call(maxbinCmd, shell=True)

    #Create contig to bin table
    bintablefile = os.path.join(binningdir_abs, 'bins_maxbin.txt')
    bintable=open(bintablefile,"a+")
    maxbindir = os.path.join(maxbinbase + '.*.fasta')
    binlist = glob.glob(maxbindir)
    for bin in binlist:
        binname = os.path.splitext(os.path.basename(bin))[0]+''
        with open(bin, 'r') as binfile:
           for line in binfile:
                if line.startswith('>'):
                    contig = line.strip()
                    contig = contig.replace(">", "")
                    bintable.write("{0}\t{1}\r\n".format(contig,binname))
    bintable.close()

    #########################
    ######## Vamb #########
    #########################

    #maxbindir = os.path.join(binningdir_abs, 'vamb')
    #if not os.path.exists(maxbindir):
        #os.makedirs(maxbindir)

    #Run vamb
    #logfile=open(logfilepath,"a+")
    #current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    #logfile.write("{0} |    Running vamb \r\n".format(current_time))
    #logfile.close()
    #maxbinCmd = 'module load anaconda3/4.4.0 vamb/20181215 && vamb -contig '+assemblypath+' -abund '+maxbindepthfile+' -out '+maxbindir+' -thread '+threads+''
    #subprocess.check_call(maxbinCmd, shell=True)

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
    bincontig_tables = ",".join(glob.glob(os.path.join(projectpath,'merged/binning', 'bins_*.txt')))
    reassemblypath = os.path.join(projectpath, 'merged', 'reassembly.fna')
    dastoolpath = os.path.join(projectpath, 'merged','binning','dastool')
    if not os.path.exists(dastoolpath):
        os.makedirs(dastoolpath)

    dastoolDependencies = 'module load gcc/5.4.0 intel/perflibs/2018 R/3.6.1 ruby/2.6.3 pullseq/1.0.2 perl/5.24.0 ncbi-blast/2.6.0+ prodigal/2.6.3 das_tool/1.1.1 diamond/0.9.24'
    dastoolCmd = ''+dastoolDependencies+' && DAS_Tool -i '+bincontig_tables+' -c '+reassemblypath+' -o '+dastoolpath+' --search_engine diamond -t '+threads+''
    subprocess.check_call(dastoolCmd, shell=True)
