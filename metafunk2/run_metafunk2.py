
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
args = parser.parse_args()

name = args.name
read1 = args.read1
read2 = args.read2
refgenpath = args.refgenpath
outpath = args.outpath
threads = args.threads

logfilepath=os.path.join(outpath,name + '.log')
logfile=open(logfilepath,"w+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
logfile.write("'{0}' | This is metafunk2 starting to run \r\n".format(current_time))
logfile.close()

#####
# Quality filtering step
#####

logfile=open(logfilepath,"a+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
logfile.write("'{0}' | This is metafunk2 starting quality filtering \r\n".format(current_time))
logfile.close()

from quality_filtering import quality_filtering
quality_filtering(read1,read2,outpath,name,threads)

#####
# Duplicate removal  step
#####

from duplicate_removal import duplicate_removal
duplicate_removal(read1,read2,outpath,name,threads)
