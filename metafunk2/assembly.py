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

def assembly(outpath,name,logfilepath,statsfilepath,threads,memory,keep,assembler):

    prevdir = "genome_mapping"
    absprevdirr = os.path.join(outpath, name + '.' + prevdir)
    #Create assembly directory if it does not exist
    assembly_dir = "assembly"
    assembly_abs = os.path.join(outpath, name + '.' + assembly_dir)
    #Create genomes directory if it does not exist
    if not os.path.exists(assembly_abs):
        os.makedirs(assembly_abs)
    #Declare input files
    read1in = os.path.join(outpath, name +  '.1.fq.gz')
    read2in = os.path.join(outpath, name +  '.2.fq.gz')

    #Input uncompressed files if compressed do not exist
    if not os.path.exists(read1in):
        read1in = os.path.join(outpath, name +  '.1.fq')
    if not os.path.exists(read2in):
        read2in = os.path.join(outpath, name +  '.2.fq')

    #Log to file
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2 has started the metagenomic assembly using {1} \r\n".format(current_time,assembler))
    logfile.close()

    #Select assembler
    if assembler == 'spades':
        assemblyCmd = 'module load tools anaconda3/2.1.0 spades/3.13.1 perl/5.20.2 && rm -rf '+assembly_abs+' && metaspades.py -1 '+read1in+' -2 '+read2in+' -t '+threads+' -m '+memory+' -k 21,29,39,59,79,99,119 --only-assembler -o '+assembly_abs+''
    if assembler == 'megahit':
        assemblyCmd = 'module load tools megahit/1.1.1 && rm -rf '+assembly_abs+' && megahit -1 '+read1in+' -2 '+read2in+' -t '+threads+' --k-list 21,29,39,59,79,99,119,141 -o '+assembly_abs+''

    #Run assembly
    subprocess.check_call(assemblyCmd, shell=True)

    #Move reads to working directory
    if assembler == 'spades':
        assembly = os.path.join(assembly_abs, 'contigs.fasta')
    if assembler == 'megahit':
        assembly = os.path.join(assembly_abs, 'final.contigs.fa')

    assemblyfinal = os.path.join(outpath, name +  '.fna')
    shutil.copy(assembly, assemblyfinal)

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
    if ( not os.path.exists(assemblyfinal) or os.stat(assemblyfinal).st_size == 0 ):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} | ERROR! Metafunk2 has stopped due to an error. Check error file \r\n".format(current_time))
        logfile.close()
        os.kill(os.getpid(), signal.SIGSTOP)
