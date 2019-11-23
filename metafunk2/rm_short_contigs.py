#!/usr/bin/env python
# -*- coding: utf-8

"""The script to assemble metagenomic reads"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip
import shutil


def rm_short_contigs(minlength,assemblyin,assemblyout):
    file = open(sys.argv[3],'a+')
    for line in open(sys.argv[2]):
        if noy line.startswith(">"): file.write(line.strip())
	    else:
            if int(line.split("_")[3])<int(sys.argv[1]): break
            else: file.write(line.strip())
    file.close()

def assembly(outpath,name,logfilepath,statsfilepath,threads,memory,keep):
    prevdir = "genome_mapping"
    absprevdirr = os.path.join(outpath, name + '.' + prevdir)
    #Create assembly directory if it does not exist
    assembly_dir = "assembly"
    assembly_abs = os.path.join(outpath, name + '.' + assembly_dir)
    #Create genomes directory if it does not exist
    if not os.path.exists(assembly_abs):
        os.makedirs(assembly_abs)
    #Declare input files
    read1in = os.path.join(outpath, name +  '.1.fq')
    read2in = os.path.join(outpath, name +  '.2.fq')

    #Run assembly
    assemblyCmd = 'metaspades.py -1 '+read1in+' -2 '+read2in+' -t '+threads+' -m '+memory+' -k 21,29,39,59,79,99,119 --only-assembler -o '+assembly_abs+''
    subprocess.check_call(assemblyCmd, shell=True)

    #Move reads to working directory
    #assembly = os.path.join(assembly_abs, 'contigs.fasta')
    #assemblyfinal = os.path.join(outpath, name +  '.fna')
    #shutil.copy(assembly, assemblyfinal)

    rm_short_contigs(1500,assembly,assemblyfinal)

    #Print error to log file if final files are not created
    if not os.path.exists(assemblyfinal):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    There was an error during the assembly. Check error file. \r\n".format(current_time,refgenname))
        logfile.close()

    #Get stats
    contigs = len([1 for line in open(assemblyfinal) if line.startswith(">")])

    #Print stats to stats file
    statsfile=open(statsfilepath,"a+")
    statsfile.write("Assembly contigs\t{0} \r\n".format(contigs))
    statsfile.close()

    #Print stats to logfile
    logfile=open(logfilepath,"a+")
    logfile.write("                     The assembly produced {0} contigs\r\n".format(contigs))
    logfile.close()

    #If keep is not selected, remove previous directory
    if not keep:
        if os.path.exists(absprevdirr):
            shutil.rmtree(absprevdirr)

    if not keep:
        if os.path.exists(assembly_abs):
            shutil.rmtree(assembly_abs)

    #Doublecheck everything is ok
    if not os.path.exists(assemblyfinal):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} | ERROR! Metafunk2 has stopped due to an error. Check error file \r\n".format(current_time))
        logfile.close()
        os.kill(os.getpid(), signal.SIGSTOP)
