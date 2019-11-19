
#Import libraries
import subprocess
import sys
import os
import argparse

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

logfile=open("log.txt","w+")
logfile.write("This is metafunk2 starting to run \r\n")

#####
# Quality filtering step
#####

from quality_filtering import quality_filtering
quality_filtering(read1,read2,outpath,name,threads)
