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
import matplotlib.pyplot as plt


def read_header(filename):
    header = {}
    with open(filename, 'r') as file:
        file.readline()  # Skip first uninform line in file ### HEADER
        for i in range(11):
            key, value = file.readline().rstrip().split(': ')
            header[key] = value
    return header


files = sorted(glob.glob('data4test/*'))
files_imps = sorted(glob.glob('data4test/c*'))

print(read_header(files[0]))
prf = np.genfromtxt(files[0], skip_header=14).T
prfimps = np.genfromtxt(files_imps[0], skip_header=14).T

seque_imps = np.hstack(prfimps[1:])
plt.close()
plt.plot(seque_imps)
plt.show()
