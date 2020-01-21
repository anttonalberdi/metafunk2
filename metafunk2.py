
# +add bin refinement
#Import libraries
import subprocess
import sys
import os
import signal
import argparse
import time
import gzip

############################
##### Argument parsing #####
############################

#Argument parsing
parser = argparse.ArgumentParser(description='Runs metafunk2 pipeline.')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-n', help="Sample name", dest="name", required=True)
required.add_argument('-1', help="Path to forward read", dest="read1", required=True)
required.add_argument('-2', help="Path to reverse read", dest="read2", required=True)
required.add_argument('-r', help="Reference genome (RG) sequence name(s) and path(s). RG1name=RG1path,RG2name=RG2path, etc.", dest="refgen", required=True)
required.add_argument('-o', help="Output path", dest="outpath", required=True)
optional.add_argument('-t', help="Number of threads (Default = 1)", dest="threads")
optional.add_argument('-m', help="RAM memory limit (Default = 50GB)", dest="memory")
optional.add_argument('-s', help="Skip steps (Default = none)", dest="skipsteps", type=str)
optional.add_argument('-i', help="Include steps (Default = all)", dest="includesteps", type=str)
optional.add_argument('-k', help="Keep intermediate files", dest="keep", action='store_true')
optional.add_argument('-a', help="Assembler software, either 'spades' or 'megahit'", dest="assembler",  type=str)
args = parser.parse_args()

name = args.name
read1 = args.read1
read2 = args.read2
refgen = args.refgen
outpath = args.outpath
threads = args.threads
memory = args.memory
keep = args.keep
assembler = args.assembler

#Prepare reference genomes
refgenlist = [l.split('=') for l in refgen.split(',') if l]
refgencount = len(refgenlist)

#Prepare threads
if args.threads is None:
    threads = 1

#Prepare memory
if args.memory is None:
    memory = 50

#Prepare skipsteps
if args.skipsteps is None:
    skipsteps = 0,
else:
    skipsteps = tuple([int(x) for x in args.skipstepsx.split(',')])
    if isinstance(skipsteps,int):
        skipsteps = (skipsteps,)

#Prepare includesteps
if args.includesteps is None:
    includesteps = 1,2,3,4,5,6,7,8,9
else:
    includesteps = tuple([int(x) for x in args.includesteps.split(',')])
    if isinstance(includesteps,int):
        includesteps = (includesteps,)

#Declare assembler
if args.assembler is None:
    assembler = 'megahit'

#############################
##### Initiate log file #####
#############################

logfilepath=os.path.join(outpath,name + '.log')
logfile=open(logfilepath,"w+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
logfile.write("\r\n\r\n====== Metafunk2 v2.0.0 ======\r\n\r\nSettings:\r\n  Threads: {1}\r\n    Memory: {2}\r\nInput files:\r\n  Read1: {3}\r\n  Read2: {4}\r\nReference genomes:\r\n    Number of reference genomes: {5} \r\n \r\n".format(current_time,threads,memory,read1,read2,refgencount))
logfile.close()

###############################
##### Initiate stats file #####
###############################

statsfilepath=os.path.join(outpath,name + '.stats')
statsfile=open(statsfilepath,"w+")
statsfile.write("Statistic\tValue \r\n".format(current_time))
statsfile.close()

##############################
#### Pipeline starts here ####
##############################

logfile=open(logfilepath,"a+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
logfile.write("\r\n{0} | Metafunk2 pipeline begins:\r\n".format(current_time))
logfile.close()

#####
## 1) QUALITY FILTERING
#############################

if ( 1 in includesteps and 1 not in skipsteps ):
    from metafunk2 import quality_filtering
    quality_filtering.quality_filtering(read1,read2,outpath,name,threads,statsfilepath,logfilepath)

#####
# 2) DUPLICATE REMOVAL
#############################

#None of the following softeare performed better thatn seqkit rmdup
##Pardre:
##https://academic.oup.com/bioinformatics/article/32/10/1562/1743431
#Fulcrum (parallel)
##https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1192-5
##FastUniq (most used)
##https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052249
##G-CNV
##https://www.frontiersin.org/articles/10.3389/fbioe.2015.00028/full

if ( 2 in includesteps and 2 not in skipsteps ):
    from metafunk2 import duplicate_removal
    duplicate_removal.duplicate_removal(read1,read2,outpath,name,threads,statsfilepath,logfilepath,keep)

#####
# 3) MAPPING AGAINST REFERENCE GENOMES
#############################

if ( 3 in includesteps and 3 not in skipsteps ):

    #Copy and prepare genomes (unzip, etc.)
    from metafunk2 import genome_mapping
    genome_mapping.copy_genome(refgenlist,outpath,name,logfilepath)

    #Index genomes if necessary
    from metafunk2 import genome_mapping
    genome_mapping.index_genome(refgenlist,outpath,name,logfilepath,threads)

    #Map sequencing reads to reference genomes
    from metafunk2 import genome_mapping
    genome_mapping.genome_mapping(refgenlist,outpath,name,logfilepath,threads,statsfilepath,keep)

#####
# 3A) Consider adding option to rarefy reads to a certain seq depth
#####

###############################################################################################################
############# THE PIPELINE SHOULD STOP HERE FOR COASSEMBLY-BASED DOWNSTREAM ANALYSES ##########################
###############################################################################################################

#####
# 4) ASSEMBLE METAGENOMIC READS
#############################

if ( 4 in includesteps and 4 not in skipsteps ):
    from metafunk2 import assembly
    assembly.assembly(outpath,name,logfilepath,statsfilepath,threads,memory,keep,assembler)

#################################################################################################################
############# THE PIPELINE SHOULD STOP HERE FOR MERGED ASSEMBLY-BASED DOWNSTREAM ANALYSES #######################
#################################################################################################################

#####
# 5) MAP READS TO ASSEMBLY
#############################

if ( 5 in includesteps and 5 not in skipsteps ):
    from metafunk2 import assembly_mapping
    assembly_mapping.assembly_mapping(outpath,name,logfilepath,threads)

#####
# 6) CONTIG BINNING
#############################

if ( 6 in includesteps and 6 not in skipsteps ):
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2 has started binning \r\n".format(current_time))
    logfile.close()

    from metafunk2 import binning
    binning.binning(outpath,name,logfilepath,threads)

    from metafunk2 import binning
    binning.bin_refinement(name,outpath,threads,memory,logfilepath)


##########################
##### Close log file #####
##########################

logfilepath=os.path.join(outpath,name + '.log')
logfile=open(logfilepath,"a+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
logfile.write("{0} | Well done! Metafunk2 has finished succesfully!\r\n\r\n".format(current_time))
logfile.close()
