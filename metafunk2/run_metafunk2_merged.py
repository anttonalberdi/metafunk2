
#Import libraries
import subprocess
import sys
import os
import argparse
import time
import gzip

#Argument parsing
parser = argparse.ArgumentParser(description='Runs metafunk2 pipeline.')
parser.add_argument('-n', help="Project name", metavar="PROJECTNAME", dest="projectname", required=True)
parser.add_argument('-p', help="Project path", metavar="PROJECTPATH", dest="projectpath", required=True)
parser.add_argument('-t', help="Number of threads", metavar="THREADS", dest="threads", default=8, required=False)
parser.add_argument('-m', help="RAM memory limit", metavar="MEMORY", dest="memory", default=250, required=False)
parser.add_argument('--skipsteps', help="Skip steps", metavar="SKIPSTEPS", dest="skipsteps", nargs='+', type=int, required=False)
parser.add_argument('--includesteps', help="Include steps", metavar="INCLUDESTEPS", dest="includesteps", nargs='+', type=int, required=False)
args = parser.parse_args()

projectname = args.projectname
projectpath = args.projectpath
threads = args.threads
memory = args.memory

#Prepare memory
if args.memory is None:
    memory = 250

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

logfilepath=os.path.join(projectpath,projectname + '.log')
logfile=open(logfilepath,"w+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
logfile.write("{0} | This is Metafunk2_merged starting to run \r\n Settings:\r\n   Threads {1}\r\n  Memory {2}\r\nInput files:\r\n   Read1 {3}\r\n  Read2 {4}\r\nReference genomes:\r\n    Number of reference genomes: {5}".format(current_time,threads,memory,read1,read2,refgencount))
logfile.close()

#####
# Initiate stats file
#####

statsfilepath=os.path.join(outpath,projectname + '.stats')
statsfile=open(statsfilepath,"w+")
statsfile.write("Statistic\tValue \r\n".format(current_time))
statsfile.close()

#####
# Checking dependencies
#####

#####
# 1) Merging assemblies
#####

if ( 1 in includesteps and 1 not in skipsteps ):
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2_merged has started to merge individual assemblies \r\n".format(current_time))
    logfile.close()

    from merge_assemblies import merge_assemblies
    merge_assemblies(projectname,projectpath,threads,memory)
