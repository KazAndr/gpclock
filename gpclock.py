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
    Help on function read_header in module gpclock:

    read_header(filename):
        return header

    Discription
    ----------
    The function reads 10 rows as a header information from file "filename".
    The function gets a name of file or path to file as input data
    and return header information.

    Parameters
    ----------
    filename : str
        Input data. Name of file in current directory or path to file.

    Returns
    -------
    header : dict
            Dictionary with header information, such as name of pulsar,
            resolution, numbers of points in the observation,
            time of start of the observation and other.

    Examples
    --------
    >> head = read_header('data4test/PULSES/compPulses_010117_1133+16_00.prf')
    >> print(head)

    {'frequency': '112.084', 'utctime': '2:19:20.921382', 'rtype': 'DPP1',
    'telcode': 'bsa1', 'obscode': 'PO', 'tau': '1.2288',
    'time': '5:19:20.92138', 'psrname': '1133+16', 'N used channels': '450',
    'date': '1/1/2017', 'period': '1.18780828035'}
    """

    header = {}
    with open(filename, 'r') as file:
        file.readline()  # Skip first uninform line in file ### HEADER
        for _ in range(11):
            key, value = file.readline().rstrip().split(': ')
            header[key] = value

    return header


def cross_correlation_using_fft(x, y):
    """
    Help on function cross_correlation_using_fft in module gpclock:

    cross_correlation_using_fft(x, y):
        return cross_function

    Discription
    ----------
        Functions calculates cross correlation function of two arrays as
        the real part of the product of the Fourier transform
        of the first series and the complex conjugate
        of the Fourier transform of the second series.

    Parameters
    ----------
    x : list, array, numpy.ndarray
        Input data. Array for correlation.

    y : list, array, numpy.ndarray
        Input data. Array for correlation.

    Returns
    -------
    cross_function : numpy.hdarray
            Cross correlation function of two arrays.

    Examples
    --------
    >> a = [0, 0, 0, 1, 0, 0]
    >> b = [0, 1, 0, 0, 0, 0]
    >> cc = cross_correlation_using_fft(a, b)
    >> print(cc)

    array([ 4.42665320e-17,  1.99863648e-17, -3.40326005e-16,  9.25185854e-18,
        1.00000000e+00, -1.40260526e-16])
    """

    f1 = fft(x)
    f2 = fft(np.flipud(y))
    cc = np.real(ifft(f1 * f2))

    return fftshift(cc)


def compute_shift(x, y):
    """
    Help on function compute_shift in module gpclock:

    compute_shift(x, y):
        return shift

    Discription
    ----------
        Functions calculates time delay between two array with
        cross correlation function.

    Parameters
    ----------
    x : list, array, numpy.ndarray
        Input data. First array.

    y : list, array, numpy.ndarray
        Input data. Second array.

    Returns
    -------
    shift : numpy.int64
            Time delay between two arrays.

    Examples
    --------
    >> a = [0, 0, 0, 1, 0, 0]
    >> b = [0, 1, 0, 0, 0, 0]
    >> compute_shift(a,b)
    -2
    >> compute_shift(b,b)
    0
    >> compute_shift(b,a)
    2
    """

    assert len(x) == len(y)
    c = cross_correlation_using_fft(x, y)
    assert len(c) == len(x)
    zero_index = int(len(x) / 2) - 1
    shift = zero_index - np.argmax(c)

    return shift


def delay_resizer(array1, array2):
    """
    Help on function dalay_resizer in module gpclock:

    dalay_resizer(array1, array2):
        return delay

    Discription
    ----------
        Functions adds zero points in the end to smallest array and caltulates
        time delay between equal size arrays.

    Parameters
    ----------
    array1 : list, array, numpy.ndarray
        Input data. First array.

    array2 : list, array, numpy.ndarray
        Input data. Second array.

    Returns
    -------
    delay : numpy.int64
            Time delay between two arrays.

    Examples
    --------
    >> a = [0, 0, 0, 1, 0, 0, 0, 0]
    >> b = [0, 1, 0, 0, 0, 0]
    >> delay_resizer(a, b)
    -2
    >> delay_resizer(b, b)
    0
    >> delay_resizer(b, a)
    2
    """

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
    Help on function get_time_delay in module gpclock:

    get_time_delay(ts1, ts2, array1, array2, tay):
        return delay

    Discription
    ----------
        Function calculates time delay between two different size array,
        witch were recorded with different time start and with
        time tesolution = tay.

    Parameters
    ----------
    ts1 : int, float
        Input data. Time start for first array.

    ts2 : int, float
        Input data. Time start for second array.

    array1 : list, array, numpy.ndarray
        Input data. First array.

    array2 : list, array, numpy.ndarray
        Input data. Second array.

    tay : int, float
        Input data. Time resolution.

    Returns
    -------
    delay : numpy.int64
            Time delay between two arrays.

    Examples
    --------
    >> a = [0, 0, 0, 1, 0, 0, 0, 0]
    >> b = [0, 1, 0, 0, 0, 0]
    >> get_time_delay(0, 0, a, b, 1)
    -2
    >> get_time_delay(0, 0, b, b, 1)
    0
    >> get_time_delay(0, 0, b, a, 1)
    2
    >>get_time_delay(2, 0, b, a, 1)
    0
    """

    if ts1 == ts2:
        delay = delay_resizer(array1, array2)
        delay *= tay

    elif ts1 > ts2:
        delay = delay_resizer(array1, array2)
        delay *= tay
        delay -= ts1 - ts2

    elif ts1 < ts2:
        delay = delay_resizer(array1, array2)
        delay *= tay
        delay += ts2 - ts1

    return delay


def get_TB_sec(filename, MJD):
    """
    Help on function get_TB_sec in module gpclock:

    get_TB_sec(file, MJD):
        return TB_sec

    Discription
    ----------
        Function returns number of seconds in
        Example of tim file can be find in data4test.

    Parameters
    ----------
    filename : str
        Input data. Name of file in current directory or path to file.

    MJD : int, str
        Input data. Modified Julian Date of time stast of observation.

    Returns
    -------
    TB_sec : nump:y.int64
            Time delay between two arrays.

    Examples
    --------
    >> get_TB_sec('data4test/_tim.out', 46436)
    58982.2204971
    """

    MJD = str(MJD)
    with open(filename, 'r') as f:
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
