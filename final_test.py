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
from gpclock import read_header, get_time_delay

earth_list = sorted(glob.glob('./final_test/*earth.csv'))
moon_list = sorted(glob.glob('./final_test/*moon.csv'))

for i, _ in enumerate(earth_list):
    header = read_header(earth_list[i])
    # print(header)
    prf_earth = np.genfromtxt(earth_list[i], skip_header=14).T

    day, month, year = header['date'].split('/')
    hour, minute, second = header['time'].split(':')
    second, microsecond = second.split('.')
    time_start_earth = datetime.datetime(
        int(year),
        int(month),
        int(day),
        int(hour),
        int(minute),
        int(second),
        int(microsecond))

    header = read_header(moon_list[i])
    prf_moon = np.genfromtxt(moon_list[i], skip_header=14).T
    # print(header)

    day, month, year = header['date'].split('/')
    hour, minute, second = header['time'].split(':')
    second, microsecond = second.split('.')
    time_start_moon = datetime.datetime(
        int(year),
        int(month),
        int(day),
        int(hour),
        int(minute),
        int(second),
        int(microsecond))

    if time_start_moon > time_start_earth:
        dt = (time_start_moon - time_start_earth)
        dt = float(str(dt.seconds) + '.' + str(dt.microseconds))*1000/float(header['tau'])
    else:
        dt = (time_start_earth - time_start_moon)
        dt = float(str(dt.seconds) + '.' + str(dt.microseconds))*1000/float(header['tau'])

    print(round(dt, 0), get_time_delay(0, 0, prf_earth[1], prf_moon[1]))
