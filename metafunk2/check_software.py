#!/usr/bin/env python
# -*- coding: utf-8

"""Checks if software is installed"""

import os
import sys
import random
import argparse
from distutils.spawn import find_executable

def is_tool(name):
    return find_executable(name) is not None
