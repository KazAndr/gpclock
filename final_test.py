#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: final_test.py
Author: Kazantsev Andrey
Email: kaz.prao@bk.ru
Github: https://github.com/KazAndr
Description: final test
"""

import glob
#import datetime

import numpy as np

from gpclock import isot_time, get_isot
from gpclock import read_header, get_time_delay, get_TB_sec, read_prf

earth_list = sorted(glob.glob('./final_test/*earth.csv'))
moon_list = sorted(glob.glob('./final_test/*moon.csv'))

def save_tim(file, time_array):
    with open(file, 'w') as f:
        for idx, value in enumerate(time_array):
            fs_p, sc_p = str(value.to_mjd()).split('.')
            sec = str(round(float( '0.' + sc_p)*24*60*60, 15))
            f.write(fs_p + '\t')
            f.write(sec + '\t')
            f.write('0.0' + '\t')
            f.write(sec + '\t')
            f.write('0.0' + '\n')

    return None

time_list_earth = []
time_list_moon = []

for i, _ in enumerate(earth_list):
    header = read_header(earth_list[i])
    # print(header)
    prf_earth = np.genfromtxt(earth_list[i], skip_header=14).T

    time_start_earth = isot_time(get_isot(header))
    time_list_earth.append(time_start_earth)

for i, _ in enumerate(moon_list):
    header = read_header(moon_list[i])
    # print(header)
    prf_earth = np.genfromtxt(moon_list[i], skip_header=14).T

    time_start_moon = isot_time(get_isot(header))
    time_list_moon.append(time_start_moon)

save_tim('./earth.out', time_list_earth)
save_tim('./moon.out', time_list_moon)

for i, _ in enumerate(earth_list):
    header = read_header(earth_list[i])
    # print(header)
    prf_earth = np.genfromtxt(earth_list[i], skip_header=14).T

    time_start_earth = isot_time(get_isot(header))
    fs_p_earth = np.round(time_start_earth.to_mjd(), 6)

    header = read_header(moon_list[i])
    prf_moon = np.genfromtxt(moon_list[i], skip_header=14).T
    # print(header)

    time_start_moon = isot_time(get_isot(header))
    fs_p_moon = np.round(time_start_moon.to_mjd(), 6)

    tay = 1.2288/1000.
    ts_e = get_TB_sec('./earth.out', fs_p_earth)
    ts_m = get_TB_sec('./moon.out', fs_p_moon)

    print(get_time_delay(0, 0, prf_earth[1], prf_moon[1], 1))
    print('-'*10)
    print(get_time_delay(ts_e, ts_m, prf_earth[1], prf_moon[1], tay))
