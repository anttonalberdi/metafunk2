
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
parser.add_argument('-r', help="Path to reference genome sequence",metavar="REFGENPATH", dest="refgenpath", required=True)
parser.add_argument('-o', help="Output path", metavar="OUTPATH", dest="outpath", required=True)
parser.add_argument('-t', help="Number of threads", metavar="THREADS", dest="threads", required=True)
parser.add_argument('--skipsteps', help="Skip steps", metavar="SKIPSTEPS", dest="skipsteps", required=False)
parser.add_argument('--includesteps', help="Include steps", metavar="INCLUDESTEPS", dest="includesteps", required=False)
args = parser.parse_args()

name = args.name
read1 = args.read1
read2 = args.read2
refgenpath = args.refgenpath
outpath = args.outpath
threads = args.threads

if args.skipsteps is None:
    skipsteps = 0
else:
    skipsteps = args.skipsteps

print(skipsteps)

if args.includesteps is None:
    includesteps = 1,2,3,4,5,6,7,8,9
else:
    includesteps = args.includesteps

print(includesteps)

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

if (1 in includesteps and 1 not in skipsteps):
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | This is metafunk2 starting quality filtering \r\n".format(current_time))
    logfile.close()

    from quality_filtering import quality_filtering
    quality_filtering(read1,read2,outpath,name,threads)

#####
# Duplicate removal step
#####

logfile=open(logfilepath,"a+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
logfile.write("{0} | This is metafunk2 starting duplicate removal \r\n".format(current_time))
logfile.close()

from duplicate_removal import duplicate_removal
duplicate_removal(read1,read2,outpath,name,threads)

#####
# Mapping against reference genomes
#####

logfile=open(logfilepath,"a+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
logfile.write("{0} | This is metafunk2 starting map reads agains reference genomes \r\n".format(current_time))
logfile.close()

from genome_mapping import copy_genome
copy_genome(refgenpath)

#from genome_mapping import check_genome
#check_genome(refgenpath)

#from genome_mapping import index_genome
#index_genome(refgenpath)

#from genome_mapping import genome_mapping
#index_genome(refgenpath)
