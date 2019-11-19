#!/usr/bin/env python
# -*- coding: utf-8

"""The script to remove duplicated sequences"""

import os
import sys
import random
import argparse
import subprocess

def duplicate_removal(read1,read2,outpath,name,threads):
    #Create quality_filtering subdirectory
    prevdir = "quality_filtering"
    absprevdirr = os.path.join(outpath, prevdir)
    newdir = "duplicate_removal"
    absnewdir = os.path.join(outpath, newdir)
    if not os.path.exists(absnewdir):
        os.makedirs(absnewdir)

    #Declare input and output files
    read1in = os.path.join(absprevdirr, name +  '.1.fq')
    read2in = os.path.join(absprevdirr, name +  '.2.fq')
    read1out = os.path.join(absprevdirr, name +  '.1.fq')
    read2out = os.path.join(absprevdirr, name +  '.2.fq')

    #Run seqkit rmdup
    Rmdup1Cmd = 'cat '+read1in+' | seqkit rmdup -s -o '+read1out+''
    Rmdup2Cmd = 'cat '+read2in+' | seqkit rmdup -s -o '+read2out+''
