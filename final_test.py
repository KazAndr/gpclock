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

from jdutil import *
from gpclock import read_header, get_time_delay, get_TB_sec, read_prf

earth_list = sorted(glob.glob('./final_test/*earth.csv'))
moon_list = sorted(glob.glob('./final_test/*moon.csv'))

def save_tim(file, time_array):
    with open(file, 'w') as f:
        for idx, value in enumerate(time_array):
            fs_p, sc_p = str(value.to_mjd()).split('.')
            sec = str(round(float( '0.' + sc_p)*24*60*60), 10**8)
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

    day, month, year = header['date'].split('/')
    hour, minute, second = header['time'].split(':')
    second, microsecond = second.split('.')
    time_start_earth = datetime(
        int(year),
        int(month),
        int(day),
        int(hour),
        int(minute),
        int(second),
        int(microsecond))
    time_list_earth.append(time_start_earth)
    
for i, _ in enumerate(moon_list):
    header = read_header(moon_list[i])
    # print(header)
    prf_earth = np.genfromtxt(moon_list[i], skip_header=14).T

    day, month, year = header['date'].split('/')
    hour, minute, second = header['time'].split(':')
    second, microsecond = second.split('.')
    time_start_moon = datetime(
        int(year),
        int(month),
        int(day),
        int(hour),
        int(minute),
        int(second),
        int(microsecond))
    time_list_moon.append(time_start_moon)

save_tim('./earth.out', time_list_earth)
save_tim('./moon.out', time_list_moon)

for i, _ in enumerate(earth_list):
    header = read_header(earth_list[i])
    # print(header)
    prf_earth = np.genfromtxt(earth_list[i], skip_header=14).T

    day, month, year = header['date'].split('/')
    hour, minute, second = header['time'].split(':')
    second, microsecond = second.split('.')
    time_start_earth = datetime(
        int(year),
        int(month),
        int(day),
        int(hour),
        int(minute),
        int(second),
        int(microsecond))
    fs_p_earth, _ = str(time_start_earth.to_mjd()).split('.')
    
    header = read_header(moon_list[i])
    prf_moon = np.genfromtxt(moon_list[i], skip_header=14).T
    # print(header)

    day, month, year = header['date'].split('/')
    hour, minute, second = header['time'].split(':')
    second, microsecond = second.split('.')
    time_start_moon = datetime(
        int(year),
        int(month),
        int(day),
        int(hour),
        int(minute),
        int(second),
        int(microsecond))
    fs_p_moon, _ = str(time_start_moon.to_mjd()).split('.')
    
    tay = 1.2288/1000.
    ts_e = get_TB_sec('./earth_0.out', fs_p_earth)
    ts_m = get_TB_sec('./moon_0.out', fs_p_moon)
    
    print(get_time_delay(0, 0, prf_earth[1], prf_moon[1], 1))
    print('-'*10)
    print(get_time_delay(ts_e, ts_m, prf_earth[1], prf_moon[1], tay))