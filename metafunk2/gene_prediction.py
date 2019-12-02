#!/usr/bin/env python
# -*- coding: utf-8

"""The script to predict genes from contigs"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip
import signal

def gene_prediction(outpath,name):
    newdir = "gene_prediction"
    absnewdir = os.path.join(outpath, name + '.' + newdir)
    if not os.path.exists(absnewdir):
        os.makedirs(absnewdir)

    contigpath = os.path.join(outpath, name + '.assembly','contigs.fasta')
    gffpath = os.path.join(absnewdir, name + '.gff')
    faapath = os.path.join(absnewdir, name + '.faa')
    fnapath = os.path.join(absnewdir, name + '.fna')

    #Run gene prediction
    predictCmd = 'prodigal -p meta -q -i '+contigpath+' -f gff -o '+gffpath+' -a '+faapath+' -d '+fnapath+''
    subprocess.check_call(predictCmd, shell=True)
