#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 12:08:56 2019

@author: Kazantsev Andrey, kaz.prao@bk.ru
"""
import os
import sys
import glob

import numpy as np

def read_header(filename):
    header = {}
    with open(filename, 'r') as file:
        file.readline() # Skip first uninform line in file ### HEADER
        for i in range(11):
            key, value = file.readline().rstrip().split(': ')
            header[key] = value
    return header

files = sorted(glob.glob('data4test/*'))

print(read_header(files[0]))
