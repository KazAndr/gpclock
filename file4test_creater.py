#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: file4test_creater.py
Author: Kazantsev Andrey
Email: kaz.prao@bk.ru
Github: https://github.com/KazAndr
Description: file for creating data for testing
"""

import os
import glob
import datetime
import numpy as np

import matplotlib.pyplot as plt
from gpclock import read_header

def wright_file(head, fName, array):
    """docstring for wright_file"""


    head_file = '### HEADER' + '\n' + \
        'telcode: ' + head['telcode'] + '\n' + \
        'obscode: ' + head['obscode'] + '\n' + \
        'rtype: ' + head['rtype'] + '\n' + \
        'psrname: ' + head['psrname'] + '\n' + \
        'period: ' + head['period'] + '\n' + \
        'tau: ' + head['tau'] + '\n' + \
        'date: ' + head['date'] + '\n' + \
        'time: ' + head['time'] + '\n' + \
        'utctime: ' + head['utctime'] + '\n' + \
        'frequency: ' + head['frequency'] + '\n' + \
        'N used channels: ' + head['N used channels'] + '\n' + \
        '### COMPENSATED PROFILE FOR EACH IMPULSE' + '\n' + \
        'time' + '\t' + 'signal'

    np.savetxt(fName,
               np.c_[range(len(array)), array],
               fmt='%i\t%1.3f',
               delimiter='\t',
               newline='\n',
               header=head_file,
               comments='')

    return None


FILES_IMP = sorted(glob.glob('data4test/PULSES/*'))

for name in FILES_IMP:
    shift = np.random.randint(-1000, 3000)
    noise_coeff = 1/2.5

    header = read_header(name)
    prfimps = np.genfromtxt(name, skip_header=14).T

    seque_imps = np.hstack(prfimps[1:])

    max_i = np.argmax(seque_imps)
    part_gp = seque_imps[max_i - 2000:max_i + 4000]

    need_points = 10000
    start_noise = np.random.normal(
        np.median(part_gp),
        noise_coeff*np.std(part_gp),
        need_points)

    end_noise = np.random.normal(
        np.median(part_gp),
        noise_coeff*np.std(part_gp),
        need_points)

    seque_imps_pl_noise_1 = np.hstack([start_noise, part_gp, end_noise])

    need_points = 10000
    start_noise = np.random.normal(
        np.median(part_gp),
        noise_coeff*np.std(part_gp),
        need_points)

    end_noise = np.random.normal(
        np.median(part_gp),
        noise_coeff*np.std(part_gp),
        need_points)

    seque_imps_pl_noise_2 = np.hstack([start_noise, part_gp, end_noise])
    seque_imps_pl_noise_2 = np.roll(seque_imps_pl_noise_2, shift)
    # # !!! Две строки, три столбца.
    # # !!! Текущая ячейка - 1
    # plt.close()
    # plt.subplot(2, 1, 1)
    # plt.plot(seque_imps_pl_noise_1)
    # plt.title("Earth pulse")
    #
    # # !!! Две строки, три столбца.
    # # !!! Текущая ячейка - 2
    # plt.subplot(2, 1, 2)
    # plt.plot(seque_imps_pl_noise_2)
    # plt.title("Moon pulse " + "Shift: " + str(shift))
    #
    # plt.tight_layout()
    #
    # plt.savefig('./out_plots/'
    #             + header['date'].replace('/', '.')
    #             + '_'
    #             + header['psrname']
    #             + '.png',
    #             format='png')
    #
fName = '.' + os.sep + os.path.basename(name)[:-4] + '.csv'
wright_file(header, fName, seque_imps_pl_noise_1)
day, month, year = header['date'].split('/')
hour, minute, second = header['time'].split(':')
second, microsecond = second.split('.')
print(datetime.datetime(
    int(year),
    int(month),
    int(day),
    int(hour),
    int(minute),
    int(second),
    int(microsecond)))
