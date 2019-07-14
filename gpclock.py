#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Jul  3 12:08:56 2019

@author: Kazantsev Andrey, kaz.prao@bk.ru
"""
# import os
# import sys
import glob
import numpy as np
import matplotlib.pyplot as plt


def read_header(filename):
    """
    Description: Read first 10 lines as header and save to dictionary
    """
    header = {}
    with open(filename, 'r') as file:
        file.readline()  # Skip first uninform line in file ### HEADER
        for _ in range(11):
            key, value = file.readline().rstrip().split(': ')
            header[key] = value
    return header


def get_time_delay(array1, array2):
    """
    Description: return corfunction as result
    """
    corfunc = np.correlate(array1, array2, 'full')
    return corfunc


FILES = sorted(glob.glob('data4test/*'))
FILES_IMP = sorted(glob.glob('data4test/c*'))

print(read_header(FILES[0]))
prf = np.genfromtxt(FILES[0], skip_header=14).T
prfimps = np.genfromtxt(FILES_IMP[0], skip_header=14).T

seque_imps = np.hstack(prfimps[1:])

plt.close()
plt.plot(seque_imps)
plt.axhline(3*np.std(seque_imps) + np.median(seque_imps))
plt.show()

need_points = 10000
start_noise = np.random.normal(
    np.mean(seque_imps),
    np.std(seque_imps),
    need_points)

end_noise = np.random.normal(
    np.mean(seque_imps),
    np.std(seque_imps),
    need_points)

seque_imps_pl_noise = np.hstack([start_noise, seque_imps, end_noise])

plt.close()
plt.plot(seque_imps_pl_noise)
# plt.axhline(3*np.std(seque_imps) + np.median(seque_imps))
plt.show()

cf = get_time_delay(seque_imps_pl_noise, np.roll(seque_imps_pl_noise, 10))

plt.close()
plt.title(np.argmax(cf))
plt.plot(cf)
# plt.axhline(3*np.std(seque_imps) + np.median(seque_imps))
plt.show()
print(len(seque_imps_pl_noise))
