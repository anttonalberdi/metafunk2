#!/usr/bin/env python
# -*- coding: utf-8

"""The script to map read to reference genomes"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip
import shutil
import signal

#Copy reference genome to working directory
def copy_genome(refgenlist,outpath,name,logfilepath):

    #Create genome_mapping directory if it does not exist
    genome_mapping_dir = "genome_mapping"
    genome_mapping_dir_abs = os.path.join(outpath, name + '.' + genome_mapping_dir)
    #Create genomes directory if it does not exist
    if not os.path.exists(genome_mapping_dir_abs):
        os.makedirs(genome_mapping_dir_abs)
    genomes_dir = "genomes"
    genomes_dir_abs = os.path.join(outpath, genomes_dir)
    if not os.path.exists(genomes_dir_abs):
        os.makedirs(genomes_dir_abs)

    #Print to Log
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2 has started to map reads against reference genomes \r\n".format(current_time))
    logfile.write("{0} |    Preparing reference genomes \r\n".format(current_time))
    logfile.close()
    #

    #Copy genome file(s) to genomes directory if it is not already there
    refgencount = len(refgenlist)
    for i in range(refgencount):
        #Declare genome name
        refgenname = refgenlist[i][0]
        #Declare original genome path
        refgenoriginalpath = refgenlist[i][1]
        #Declare new genome path and copy file
        if refgenoriginalpath.endswith('.gz'):
            refgenpath = os.path.join(genomes_dir_abs, refgenname + '.fna.gz')
            refgenpathuncomp = refgenpath.replace(".gz", "")
            if ( not os.path.exists(refgenpath) and not os.path.exists(refgenpathuncomp) ):
                logfile=open(logfilepath,"a+")
                current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
                logfile.write("{0} |        Transferring {1} genome to working directory \r\n".format(current_time,refgenname))
                logfile.close()
                shutil.copy(refgenoriginalpath, refgenpath)
            else:
                logfile=open(logfilepath,"a+")
                current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
                logfile.write("{0} |        {1} genome is already in the working directory \r\n".format(current_time,refgenname))
                logfile.close()
        if ( refgenoriginalpath.endswith('.fasta') or refgenoriginalpath.endswith('.fa') or refgenoriginalpath.endswith('.fna') ):
            refgenpath = os.path.join(genomes_dir_abs, refgenname + '.fna')
            if not os.path.exists(refgenpath):
                logfile=open(logfilepath,"a+")
                current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
                logfile.write("{0} |        Transferring {1} genome to working directory \r\n".format(current_time,refgenname))
                logfile.close()
                shutil.copy(refgenoriginalpath, refgenpath)
            else:
                logfile=open(logfilepath,"a+")
                current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
                logfile.write("{0} |        {1} genome is already in the working directory \r\n".format(current_time,refgenname))
                logfile.close()

        #Manipulate reference genome
        if refgenpath.endswith('.gz'):
            refgenpathuncomp = refgenpath.replace(".gz", "")
            if not os.path.exists(refgenpathuncomp):
                logfile=open(logfilepath,"a+")
                current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
                logfile.write("{0} |        Decompressing {1} genome \r\n".format(current_time,refgenname))
                logfile.close()
                DecompCmd = 'module load tools pigz/2.3.4 && pigz -d '+refgenpath+''
                subprocess.check_call(DecompCmd, shell=True)

#Index reference genomes if not already indexed or being indexed by another job
def index_genome(refgenlist,outpath,name,logfilepath,threads):
    refgencount = len(refgenlist)
    for i in range(refgencount):
        #Declare genome name and path
        refgenname = refgenlist[i][0]
        refgenpath = os.path.join(outpath,'genomes', refgenname + '.fna')
        refgenflag = os.path.join(outpath,'genomes', refgenname + '.indexing')
        refgenfai = os.path.join(outpath,'genomes', refgenname + '.fna.fai')
        if ( not os.path.exists(refgenfai) and not os.path.exists(refgenflag) ):
            #Create indexing flag file to let other jobs know this genome is being indexed
            indexingflagfile=open(refgenflag,"w+")
            indexingflagfile.write("Genome {0} is being indexed".format(refgenname))
            indexingflagfile.close()
            #Print to log file
            logfile=open(logfilepath,"a+")
            current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
            logfile.write("{0} |        Indexing {1} genome \r\n".format(current_time,refgenname))
            logfile.close()
            #Index genomes
            logfile=open(logfilepath,"a+")
            current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
            logfile.write("{0} |            Samtools index \r\n".format(current_time,refgenname))
            logfile.close()
            samtoolsindexCmd = 'module load tools samtools/1.9 && samtools faidx '+refgenpath+''
            subprocess.check_call(samtoolsindexCmd, shell=True)

            logfile=open(logfilepath,"a+")
            current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
            logfile.write("{0} |            Bwa index \r\n".format(current_time,refgenname))
            logfile.close()
            bwaindexCmd = 'module load tools bwa/0.7.15 && bwa index '+refgenpath+''
            subprocess.check_call(bwaindexCmd, shell=True)

            #Hisat2 is too slow compared to bwa
            #logfile=open(logfilepath,"a+")
            #current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
            #logfile.write("{0} |        Hisat index \r\n".format(current_time,refgenname))
            #logfile.close()
            #hisatindexCmd = 'module unload gcc/5.1.0 && module load tools anaconda2/4.4.0 hisat2/2.1.0 && hisat2-build '+refgenpath+' '+refgenpath+' -p '+threads+''
            #subprocess.check_call(hisatindexCmd, shell=True)

            #Remove indexing flag when indexing is done
            os.remove(refgenflag)
        else:
            logfile=open(logfilepath,"a+")
            current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
            logfile.write("{0} |        {1} genome is already indexed \r\n".format(current_time,refgenname))
            logfile.close()


###
# add waiting https://blog.miguelgrinberg.com/post/how-to-make-python-wait
###

def genome_mapping(refgenlist,outpath,name,logfilepath,threads,statsfilepath,keep):
    #Declare source directory
    prevdir = "duplicate_removal"
    absprevdirr = os.path.join(outpath, name + '.' + prevdir)

    #Iterate across reference genomes
    refgencount = len(refgenlist)
    for i in range(refgencount):
        #Declare input files
        #If first genome, get files from previous step directory
        if i==0:
            read1in = os.path.join(absprevdirr, name +  '.1.fq')
            read2in = os.path.join(absprevdirr, name +  '.2.fq')
        #If not first genome, get files outputed from previous genome mapping
        else:
            read1in = os.path.join(outpath,name + '.genome_mapping', name +  '.1.fq')
            read2in = os.path.join(outpath,name + '.genome_mapping', name +  '.2.fq')

        #Declare genome name and file paths
        refgenname = refgenlist[i][0]
        refgenpath = os.path.join(outpath,'genomes', refgenname + '.fna')
        refgenflag = os.path.join(outpath,'genomes', refgenname + '.indexing')
        bampath_all = os.path.join(outpath, name + '.genome_mapping', name + '.' + refgenname + '.bam')
        bampath_host = os.path.join(outpath,name + '.genome_mapping', name + '.mappedto.' + refgenname + '.bam')
        bampath_mg = os.path.join(outpath,name + '.genome_mapping', name + '.mg.bam')
        read1out = os.path.join(outpath,name + '.genome_mapping', name +  '.1.fq')
        read2out = os.path.join(outpath,name + '.genome_mapping', name +  '.2.fq')

        #Wait until indexing flag disappears (meaning another job is indexing the genome)
        #The script will check the presence of the indexing file every 5 minutes and proceed when  disappears
        secs = 0
        while os.path.exists(refgenflag):
            secs = secs + 300
            mins = secs / 60
            #Print to log
            logfile=open(logfilepath,"a+")
            current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
            logfile.write("{0} |    Waiting {1} minutes until genome {2} is indexed \r\n".format(current_time,mins,refgenname))
            logfile.close()
            time.sleep(secs)

        #Declare mapping commands
        #hisat version mapCmd = 'module load tools module load tools anaconda2/4.4.0 hisat2/2.1.0 samtools/1.9 && hisat2 -x '+refgenpath+' -1 '+read1in+' -2 '+read2in+' -p '+threads+' | samtools view -T '+refgenpath+' -b - > '+bampath_all+''
        mapCmd = 'module load tools samtools/1.9 bwa/0.7.15 && bwa mem -t '+threads+' -R "@RG\tID:ProjectName\tCN:AuthorName\tDS:Mappingt\tPL:Illumina1.9\tSM:Sample" '+refgenpath+' '+read1in+' '+read2in+' | samtools view -T '+refgenpath+' -b - > '+bampath_all+''
        hostmapCmd = 'module load tools samtools/1.9 && samtools view -T '+refgenpath+' -b -F12 '+bampath_all+' > '+bampath_host+''
        mgmapCmd = 'module load tools samtools/1.9 && samtools view -T '+refgenpath+' -b -f12 '+bampath_all+' > '+bampath_mg+''
        mgfqCmd = 'module load tools samtools/1.9 && samtools fastq -1 '+read1out+' -2 '+read2out+' '+bampath_mg+''

        #Mapping to genome
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Mapping reads to {1} genome \r\n".format(current_time,refgenname))
        logfile.close()
        subprocess.check_call(mapCmd, shell=True)

        #Mapping to genome using Hisat2 - too slow
        #logfile=open(logfilepath,"a+")
        #current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        #logfile.write("{0} |    Mapping reads to {1} genome using hisat \r\n".format(current_time,refgenname))
        #logfile.close()
        #subprocess.check_call(mapCmd2, shell=True)

        #Extracting mapped (genomic) reads
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Extracting reads mapped to {1} genome (genomic reads) \r\n".format(current_time,refgenname))
        logfile.close()
        subprocess.check_call(hostmapCmd, shell=True)

        #Extracting unmapped (metagenomic) reads
        #Extracting mapped (genomic) reads
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} |    Extracting reads not mapped to {1} genome (metagenomic reads) and transforming them to fastq \r\n".format(current_time,refgenname))
        logfile.close()
        #Extract unmapped
        subprocess.check_call(mgmapCmd, shell=True)
        #Convert unmapped to fastq
        subprocess.check_call(mgfqCmd, shell=True)

        #Get stats
        reads = 0
        bases = 0
        with open(read1out, 'rb') as read:
            for id in read:
                seq = next(read)
                reads += 1
                bases += len(seq.strip())*2
                next(read)
                next(read)
        #Print stats to statsfile
        statsfile=open(statsfilepath,"a+")
        statsfile.write("Reads after mapping to {0}\t{1} ({2} bases)\r\n".format(refgenname,reads,bases))
        statsfile.close()
        #Print stats to logfile
        logfile=open(logfilepath,"a+")
        logfile.write("                         {1} reads ({2} bases) were kept after mapping to {0} genome\r\n".format(refgenname,reads,bases))
        logfile.close()

    ###
    # Per-genome iteration ends here
    ###

    #Compress and move final read files to parent folder
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Compressing metagenomic read files\r\n".format(current_time,refgenname))
    logfile.close()
    read1final = os.path.join(outpath, name +  '.1.fq.gz')
    read2final = os.path.join(outpath, name +  '.2.fq.gz')
    read1Cmd = 'module load tools pigz/2.3.4 && pigz -p '+threads+' -c '+read1out+' > '+read1final+''
    subprocess.check_call(read1Cmd, shell=True)
    read2Cmd = 'module load tools pigz/2.3.4 && pigz -p '+threads+' -c '+read2out+' > '+read2final+''
    subprocess.check_call(read2Cmd, shell=True)

    #If keep is not selected, remove previous directory
    if not keep:
        if os.path.exists(absprevdirr):
            shutil.rmtree(absprevdirr)

    #Doublecheck if everything is ok
    if ( os.stat(read1final).st_size == 0 or  os.stat(read2final).st_size == 0):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("{0} | ERROR! Metafunk2 has stopped due to an error. Check error file \r\n".format(current_time))
        logfile.close()
        os.kill(os.getpid(), signal.SIGSTOP)
