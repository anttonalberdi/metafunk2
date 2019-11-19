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
    subdir = "duplicate_removal"
    absdir = os.path.join(outpath, subdir)
    if not os.path.exists(absdir):
        os.makedirs(absdir)
