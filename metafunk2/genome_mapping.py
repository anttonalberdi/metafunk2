#!/usr/bin/env python
# -*- coding: utf-8

"""The script to map read to reference genomes"""

import os
import sys
import random
import argparse
import subprocess

def check_genome(refgenomepath):
    if refgenpath.endswith('.gz'):
        logfile=open(logfilepath,"a+")
        current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
        logfile.write("     {0} | Decompressing reference genome {1} \r\n".format(current_time,refgenomepath))
        logfile.close()
        DecompCmd = 'pigz '+refgenpath+''
        subprocess.check_call(DecompCmd, shell=True)
        refgenomepath = refgenomepath.replace(".gz", "")

def index_genome(refgenomepath):


def genome_mapping(read1,read2,refgenomepath):
