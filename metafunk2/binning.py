### BINNING UPDATE + bin refinement

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
#added glob
import glob
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

    #Create contig to bin table
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Writing bin table for metabat \r\n".format(current_time))
    logfile.close()
    bintablefile = os.path.join(absnewdir, 'bins_metabat.txt')
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

    #Create contig to bin table
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Writing bin table for maxbin \r\n".format(current_time))
    logfile.close()
    bintablefile = os.path.join(absnewdir, 'bins_maxbin.txt')
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



###################### ADD BIN TABLES TO EACH ONE, MAKE THAT IT IS SAVED if __name__ == '__main__':
## name.binning, not inside name.binning/metabat or maxbin!!!
# The input for das_tool has to be a bin table, not paths!



def bin_refinement(name,outpath,threads,memory,logfilepath):
    #              (projectname,projectpath,threads,memory,logfilepath)
#/home/projects/ku-cbd/people/nurher/chick_metafunk2_test/CA16_13F1b.binning/metabat
    print("Start bin refinement") #added

    bincontig_tables = ",".join(glob.glob(os.path.join(outpath,name,'binning','bins_*.txt'))) #CHANGED
    #Input
    assemblypath = os.path.join(outpath, name + '.fna') #CHANGED
    dastoolpath = os.path.join(outpath,name,'dastool') #CHANGED
    if not os.path.exists(dastoolpath):
        os.makedirs(dastoolpath)
    dastoolbase = os.path.join(dastoolpath, 'dastool')

    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Refinning bins using DAS_Tool \r\n".format(current_time))
    logfile.close()

    #Refinement using DAS_Tool
    dastooldb = '/home/projects/ku-cbd/people/antalb/databases/dastool_db'
    dastoolDependencies = 'module load tools gcc/5.4.0 intel/perflibs/2018 R/3.6.1 ruby/2.6.3 pullseq/1.0.2 perl/5.24.0 ncbi-blast/2.6.0+ prodigal/2.6.3 das_tool/1.1.1 diamond/0.9.24 usearch/11.0.667'
    dastoolCmd = ''+dastoolDependencies+' && DAS_Tool -i '+bincontig_table+' -c '+assemblypath+' -o '+dastoolbase+' -l maxbin,metabat --search_engine diamond -t '+threads+' --db_directory '+dastooldb+' --write_bins 1'
    subprocess.check_call(dastoolCmd, shell=True)

    #Refinement using Binning_refiner (problems with R dependencies)
    #module unload gcc gcc/5.1.0 && module load anaconda3/4.0.0 && Binning_refiner -i metafunk2_test2/merged/binning/refiner/ -p refined -plot

    #Move definitive bins to binning directory
    binsource = os.path.join(dastoolpath)   # dastoolpath,'dastool_DASTool_bins')
    bindestination = os.path.join(outpath,name,'binning') #CHANGED
    binfiles = glob.glob(os.path.join(binsource,'*.fa'))
    for b in binfiles:
        shutil.move(b, bindestination)
