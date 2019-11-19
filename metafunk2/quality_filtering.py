#!/usr/bin/env python
# -*- coding: utf-8

"""The script to quality-filter the sequencing reads"""

import os
import sys
import random
import argparse

def quality_filtering(read1,read2,outpath,name,threads):
    myCmd = 'AdapterRemoval --file1 '+read1+' --file2 '+read2+' --basename '+outpath+'/quality_filtering/'+name+' --minquality 30 --trimqualities --trimns --maxns 5 --threads '+threads+''
    os.system(myCmd)
