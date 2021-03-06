
#Import libraries
import subprocess
import sys
import os
import argparse
import time
import gzip

############################
##### Argument parsing #####
############################

parser = argparse.ArgumentParser(description='Runs metafunk2 merged pipeline.')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-n', help="Project name", metavar="PROJECTNAME", dest="projectname", required=True)
required.add_argument('-p', help="Project path", metavar="PROJECTPATH", dest="projectpath", required=True)
optional.add_argument('-t', help="Number of threads", metavar="THREADS", dest="threads", default=8, required=False)
optional.add_argument('-m', help="RAM memory limit", metavar="MEMORY", dest="memory", default=250, required=False)
optional.add_argument('-s', help="Skip steps", dest="skipsteps", type=str)
optional.add_argument('-i', help="Include steps", dest="includesteps", type=str)
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

#############################
##### Initiate log file #####
#############################

logfilepath=os.path.join(projectpath,projectname + '.log')
logfile=open(logfilepath,"w+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
logfile.write("{0} | This is Metafunk2_merged starting to run \r\n".format(current_time))
logfile.close()

###############################
##### Initiate stats file #####
###############################

statsfilepath=os.path.join(projectpath,projectname + '.stats')
statsfile=open(statsfilepath,"w+")
statsfile.write("Statistic\tValue \r\n".format(current_time))
statsfile.close()

###############################################################################################################
############# THE PIPELINE REQUIRES METAFUNK2 PIPELINE STEPS 1,2,3,4 TO BE PREVIOUSLY PERFORMED ##########################
###############################################################################################################

#####
# 1) Merging assemblies
#####

if ( 1 in includesteps and 1 not in skipsteps ):
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2_merged is merging assemblies \r\n".format(current_time))
    logfile.close()

    from metafunk2 import merge_assemblies
    merge_assemblies.merge_assemblies(projectname,projectpath,threads,memory,logfilepath)

#####
# 2) Merged assembly mapping
#####

if ( 2 in includesteps and 2 not in skipsteps ):
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2_merged is mapping reads of different samples to the reassembly (merged assemblies) \r\n".format(current_time))
    logfile.close()

    from metafunk2 import reassembly_mapping
    reassembly_mapping.reassembly_indexing(projectname,projectpath,threads,memory,logfilepath)

    from metafunk2 import reassembly_mapping
    reassembly_mapping.reassembly_mapping(projectname,projectpath,threads,memory,logfilepath)

#####
# 3) Merged assembly binning
#####

if ( 3 in includesteps and 3 not in skipsteps ):
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2_merged is binning \r\n".format(current_time))
    logfile.close()

    from metafunk2 import binning_merged
    binning_merged.binning_merged(projectname,projectpath,threads,memory,logfilepath)

    from metafunk2 import binning_merged
    binning_merged.bin_refinement(projectname,projectpath,threads,memory,logfilepath)

#####
# 4) Bin quality check
#####

if ( 4 in includesteps and 4 not in skipsteps ):
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2_merged is assessing quality of bins \r\n".format(current_time))
    logfile.close()

    from metafunk2 import bin_qc
    bin_qc.bin_qc(projectpath,threads,memory,logfilepath)

#####
# 5) Bin annotation
#####

if ( 5 in includesteps and 5 not in skipsteps ):
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} | Metafunk2_merged is annotating the bins \r\n".format(current_time))
    logfile.close()

    from metafunk2 import bin_annotation
    bin_annotation.bin_annotation(projectpath,threads,memory,logfilepath)


#####
# Close log file
#####

logfilepath=os.path.join(projectpath,projectname + '.log')
logfile=open(logfilepath,"a+")
current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
logfile.write("{0} | Well done! Metafunk2_merged has finished succesfully!\r\n\r\n".format(current_time))
logfile.close()
