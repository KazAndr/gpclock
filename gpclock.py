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

def delay_resizer(array1, array2):
    """docstring for delay_resizer"""
    if len(array1) == len(array2):
        cf = np.correlate(array1, array2, 'full')
        delay = ((len(cf)-1)/2)-np.argmax(cf)

    elif len(array1) > len(array2):
        array2 = np.append(array2, np.zeros(len(array1) - len(array2)))
        cf = np.correlate(array1, array2, 'full')
        # delay = (l_frame, r_frame)
        delay = ((len(cf)-1)/2)-np.argmax(cf)

    elif len(array1) < len(array2):
        array1 = np.append(array1, np.zeros(len(array2) - len(array1)))
        cf = np.correlate(array2, array1, 'full')
        delay = -(((len(cf)-1)/2)-np.argmax(cf))
    return delay


def get_time_delay(ts1, ts2, array1, array2):

    """
    Description: return corfunction as result
    """
    if ts1 == ts2:
        delay = delay_resizer(array1, array2)

    elif ts1 > ts2:
        array1 = np.append(np.zeros(ts1 - ts2), array1)
        delay = delay_resizer(array1, array2)

    elif ts1 < ts2:
        array2 = np.append(np.zeros(ts2 - ts1), array2)
        delay = delay_resizer(array1, array2)

    return delay


FILES = sorted(glob.glob('data4test/AP/*'))
FILES_IMP = sorted(glob.glob('data4test/PULSES/*'))

# print(read_header(FILES[0]))
prf = np.genfromtxt(FILES[0], skip_header=14).T
prfimps = np.genfromtxt(FILES_IMP[0], skip_header=14).T

seque_imps = np.hstack(prfimps[1:])

# plt.close()
# plt.plot(seque_imps)
# plt.axhline(3*np.std(seque_imps) + np.median(seque_imps))
# plt.show()

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


# plt.close()
# plt.title(np.argmax(cf))
# plt.plot(cf)
# # plt.axhline(3*np.std(seque_imps) + np.median(seque_imps))
# plt.show()
ts1, ts2 = 100, 10
for i in range(-2, 2):
    moon_pulse = seque_imps_pl_noise[ts1:-800]
    swich_moon_pulse = np.roll(seque_imps_pl_noise, i)[ts2:-122]
    # cf = get_time_delay(swich_moon_pulse, moon_pulse)
    cf = get_time_delay(ts1, ts2, moon_pulse, swich_moon_pulse)
    print(i, cf)
# plt.close()
# plt.plot(moon_pulse)
# plt.plot(swich_moon_pulse)
# plt.xlim(7500, 17500)
# plt.show()
