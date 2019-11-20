
#Import libraries
import subprocess
import sys
import os
import argparse
import time

#Argument parsing
parser = argparse.ArgumentParser(description='Runs metafunk2 pipeline.')
parser.add_argument('-n', help="Sample name",metavar="NAME", dest="name", required=True)
parser.add_argument('-1', help="KEGG database path",metavar="READ1", dest="read1", required=True)
parser.add_argument('-2', help="KEGG database path",metavar="READ2", dest="read2", required=True)
parser.add_argument('-r', help="Reference genome sequences",metavar="REFGEN", dest="refgen", required=True)
parser.add_argument('-o', help="Output path", metavar="OUTPATH", dest="outpath", required=True)
parser.add_argument('-t', help="Number of threads", metavar="THREADS", dest="threads", required=True)
parser.add_argument('--skipsteps', help="Skip steps", metavar="SKIPSTEPS", dest="skipsteps", nargs='+', type=int, required=False)
parser.add_argument('--includesteps', help="Include steps", metavar="INCLUDESTEPS", dest="includesteps", nargs='+', type=int, required=False)
args = parser.parse_args()

name = args.name
read1 = args.read1
read2 = args.read2
outpath = args.outpath
threads = args.threads

#Prepare reference genomes
refgen = args.refgen
refgenlist = [l.split('=') for l in refgen.split(',') if l]
refgencount = len(refgenlist)

print(refgen)
print(refgenlist)
print(refgencount)

#Prepare skipsteps
if args.skipsteps is None:
    skipsteps = 0,
else:
    skipsteps = tuple(args.skipsteps)
    if isinstance(skipsteps,int):
        skipsteps = (skipsteps,)

#Prepare includesteps
if args.includesteps is None:
    includesteps = 1,2,3,4,5,6,7,8,9
else:
    includesteps = tuple(args.includesteps)
    if isinstance(includesteps,int):
        includesteps = (includesteps,)

#####
# Initiate log file
#####

logfilepath=os.path.join(outpath,name + '.log')
logfile=open(logfilepath,"w+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
logfile.write("{0} | This is metafunk2 starting to run \r\n".format(current_time))
logfile.close()

#####
# Checking dependencies
#####

#####
# Quality filtering step
#####

if ( 1 in includesteps and 1 not in skipsteps ):
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2 has started quality filtering the reads \r\n".format(current_time))
    logfile.close()

    from quality_filtering import quality_filtering
    quality_filtering(read1,read2,outpath,name,threads)

#####
# Duplicate removal step
#####

if ( 2 in includesteps and 2 not in skipsteps ):
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2 has started to remove duplicated reads (clonality) \r\n".format(current_time))
    logfile.close()

    from duplicate_removal import duplicate_removal
    duplicate_removal(read1,read2,outpath,name,threads)

#####
# Mapping against reference genomes
#####

if ( 3 in includesteps and 3 not in skipsteps ):
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2 has started to map reads against reference genomes \r\n".format(current_time))
    logfile.close()

    from genome_mapping import copy_genome
    copy_genome(refgenlist,outpath,name,logfilepath)

    #from genome_mapping import edit_genome
    #edit_genome(refgenpath,outpath,name,logfilepath)

#from genome_mapping import index_genome
#index_genome(refgenpath)

#from genome_mapping import genome_mapping
#index_genome(refgenpath)
