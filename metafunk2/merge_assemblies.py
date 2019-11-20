#!/usr/bin/env python
# -*- coding: utf-8

"""The script for merging assemblies"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip

def merge_assemblies(projectname,projectpath,threads,memory):
    #Create merged directory if it does not exist
    merged_dir = "merged"
    merged_abs = os.path.join(projectpath, merged_dir)
    if not os.path.exists(merged_abs):
        os.makedirs(merged_abs)

    assembliespath = os.path.join(projectpath,'*.assembly', 'contigs.fasta')
    mergedassembliespath = os.path.join(merged_abs, 'assemblies.fna')

    #Concatenate assemblies
    concCmd = 'cat '+assembliespath+' > '+mergedassembliespath+''
    subprocess.check_call(concCmd, shell=True)
