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

from numpy.fft import fft, ifft, fft2, ifft2, fftshift

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

def cross_correlation_using_fft(x, y):
    f1 = fft(x)
    f2 = fft(np.flipud(y))
    cc = np.real(ifft(f1 * f2))
    return fftshift(cc)

# shift &lt; 0 means that y starts 'shift' time steps before x # shift &gt; 0 means that y starts 'shift' time steps after x
def compute_shift(x, y):
    assert len(x) == len(y)
    c = cross_correlation_using_fft(x, y)
    assert len(c) == len(x)
    zero_index = int(len(x) / 2) - 1
    shift = zero_index - np.argmax(c)
    return shift

def delay_resizer(array1, array2):
    """docstring for delay_resizer"""
    if len(array1) == len(array2):
        delay = compute_shift(array1, array2)

    elif len(array1) > len(array2):
        array2 = np.append(array2, np.zeros(len(array1) - len(array2)))
        delay = compute_shift(array1, array2)

    elif len(array1) < len(array2):
        array1 = np.append(array1, np.zeros(len(array2) - len(array1)))
        delay = compute_shift(array1, array2)

    return delay

def get_time_delay(ts1, ts2, array1, array2, tay):
    """
    Description: return corfunction as result
    """
    if ts1 == ts2:
        delay = delay_resizer(array1, array2)

    elif ts1 > ts2:
        delay = delay_resizer(array1, array2)
        delay -= ts1 - ts2

    elif ts1 < ts2:
        delay = delay_resizer(array1, array2)
        delay += ts2 - ts1

    return delay*tay

def get_TB_sec(file, MJD):
    MJD = str(MJD)
    with open(file, 'r') as f:
        lines = f.readlines()

    tim = [line.split() for line in lines if len(line.split()) == 7]
    tim_array = np.asarray(tim).T
    idx = np.where(tim_array[2] == MJD)

    return float(tim_array[5][idx].item())

if __name__ == "__main__":
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
    ts1, ts2 = 0, 30
    for i in range(400, 10000):
        moon_pulse = seque_imps_pl_noise[ts1:-8]
        swich_moon_pulse = np.roll(seque_imps_pl_noise, i)[ts2:-1225]
        # cf = get_time_delay(swich_moon_pulse, moon_pulse)
        cf = get_time_delay(ts1, ts2, moon_pulse, swich_moon_pulse)
        print(i, cf)
    # plt.close()
    # plt.plot(moon_pulse)
    # plt.plot(swich_moon_pulse)
    # plt.xlim(7500, 17500)
    # plt.show()
