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

from gpclock import isot_time, get_isot, read_header, get_time_delay

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

shift_list = []
cut_1_list = []
cut_2_list = []
dt_start_list = []
results_shift = []

for name in FILES_IMP:
    shift = np.random.randint(-1000, 3000)
    cut_1 = np.random.randint(-1000, 0)
    cut_2 = np.random.randint(-1000, 0)
    dt_start = np.random.randint(0, 1000)

    shift_list.append(shift)
    cut_1_list.append(cut_1)
    cut_2_list.append(cut_2)
    dt_start_list.append(dt_start)


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
    seque_imps_pl_noise_1 = seque_imps_pl_noise_1[:cut_1]

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
    seque_imps_pl_noise_2 = seque_imps_pl_noise_2[:cut_2]



    # !!! Две строки, три столбца.
    # !!! Текущая ячейка - 1
    plt.close()
    plt.subplot(2, 1, 1)
    plt.plot(seque_imps_pl_noise_1)
    plt.title("Earth pulse " + "Cut: " + str(cut_1) + " points")
    plt.ylabel("Intensity, arb.u.")
    plt.xlabel("Number of points, dt = " + header['tau'])
    # !!! Две строки, три столбца.
    # !!! Текущая ячейка - 2
    plt.subplot(2, 1, 2)
    plt.plot(seque_imps_pl_noise_2)
    plt.title("Moon pulse " + "Shift: " + str(shift) + 'points, Cut: ' + str(cut_2) + " points")
    plt.ylabel("Intensity, arb.u.")
    plt.xlabel("Number of points, dt = " + header['tau'])

    plt.tight_layout()

    plt.savefig('./out_plots/'
                + header['date'].replace('/', '.')
                + '_'
                + header['psrname']
                + '.png',
                format='png', dpi=150)


    time_start = isot_time(get_isot(header))

    if dt_start > 0:
        seque_imps_pl_noise_2 = seque_imps_pl_noise_2[abs(dt_start):]
    else:
        seque_imps_pl_noise_1 = seque_imps_pl_noise_1[abs(dt_start):]

    results_shift.append(get_time_delay(0,0,
                                        seque_imps_pl_noise_1,
                                        seque_imps_pl_noise_2, 1))

    fName = './final_test' + os.sep + os.path.basename(name)[:-4] + '_earth.csv'
    wright_file(header, fName, seque_imps_pl_noise_1)


    dt_s, dt_ms = str(round(dt_start * float(header['tau'])/1000, 6)).split('.')

    shifted_time = str(
        np.datetime64(time_start.value)
        + np.timedelta64(dt_s, 's')
        + np.timedelta64(dt_ms + '00', 'us')
    )

    shifted_time = isot_time(shifted_time)

    header['date'] = (str(shifted_time.day) + '/'
                      + str(shifted_time.month) + '/'
                      + str(shifted_time.year))

    header['utctime'] = (str(shifted_time.hour) + ':'
                      + str(shifted_time.minutes) + ':'
                      + str(shifted_time.seconds))

    fName = './final_test' + os.sep + os.path.basename(name)[:-4] + '_moon.csv'
    wright_file(header, fName, seque_imps_pl_noise_2)

fName = '.' + os.sep + 'shift.csv'
head_file = '#SHIFT   CUT1   CUT2   DT_START   RES_SHIFT'
np.savetxt(fName,
           np.c_[shift_list, cut_1_list, cut_2_list, dt_start_list, results_shift],
           fmt='%i\t%i\t%i\t%i\t%i',
           delimiter='\t',
           newline='\n',
           header=head_file,
           comments='')
