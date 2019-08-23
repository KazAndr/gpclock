#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Jul  3 12:08:56 2019

@author: Kazantsev Andrey, kaz.prao@bk.ru
"""

import numpy as np

from numpy.fft import fft, ifft, fftshift

class isot_time(object):
    """
    Descriotion on class isot_time(object) in module gpclock:

    Discription
    ----------
    The class allows to effectively work with time and date in frame the module.
    """
    def __init__(self, time):
        """
        Help on method __init__ in class isot_time in module gpclock:

        Discription
        ----------
        The method initializes the class isot_time.

        Parameters
        ----------
        time : str
            Input data. Time for initialization in isot format
            (YYYY-MM-DDThh:mm:ss.ssssss).

        Attributes:
        -------
        value : str
            Full value of input parametr 'time'.

        year : str
            Value of year in 'time'.

        month : str
            Value of month in 'time'.

        day : str
            Value of day in 'time'.

        hour : str
            Value of hour in 'time'.

        minutes : str
            Value of minutes in 'time'.

        seconds : str
            Value of seconds in 'time'.

        Examples
        --------
        >> t = isot_time('2002-02-03T13:56:03.1722583')
        >> t.value
        '2002-02-03T13:56:03.1722583'
        >> t.year
        '2002'
        >> t.minutes
        '56'
        >> t.seconds
        '03.1722583'
        """

        self.value = time
        ymd, hms = time.split('T')
        self.year, self.month, self.day = ymd.split('-')
        self.hour, self.minutes, self.seconds = hms.split(':')

    def to_mjd(self):
        """
        Help on method to_mjd in class isot_time in module gpclock:

        Discription
        ----------
            The method converts 'time' to MJD(Modified Julian Date) format
            by algoritm from (https://en.wikipedia.org/wiki/Julian_day#Variants)
            and returns value of MJD.

        Parameters
        ----------
        The method uses attributes of class(see __init__).

        Returns
        -------
        MJD : numpy.float64
                Value of MJD.

        Examples
        --------
        >> t = isot_time('2002-02-03T13:56:03.1722583')
        >> t.to_mjd()
        52308.580592271406
        """
        year = np.int(self.year)
        month = np.int(self.month)
        day = np.int(self.day)
        hh = np.int(self.hour)
        mm = np.int(self.minutes)
        ss = np.float64(self.seconds)

        a=int((14-month)/12)
        y=year+4800-a
        m=month+12*a-3
        JDN=day+int((153*m+2)/5)+int(365*y)+int(y/4)-int(y/100)+int(y/400)-32045

        JD=JDN+((hh-12)/24)+(mm/1440)+(ss/86400)

        return JD - 2400000.5


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


def read_prf(filename):
    """
    Help on function read_prf in module gpclock:

    read_prf(filename):
        return header, time_point, intens_point

    Discription
    ----------
    The function returns header information, time array and flux array form
    file of observations.

    Parameters
    ----------
    filename : str
        Input data. Name of file in current directory or path to file.

    Returns
    -------
    header : dict
        See documentation for read_header function.

    time_point : numpy.ndarray
        Array of time points from file of observations.

    intens_point : numpy.ndarray
        Array of values of intensity from file of observations.
    Examples
    --------
    >> head, t_point, i_point = read_prf('data4test/AP/010117_1133+16_00.prf')
    >> print(head)

    {'frequency': '112.084', 'utctime': '2:19:20.921382', 'rtype': 'DPP1',
    'telcode': 'bsa1', 'obscode': 'PO', 'tau': '1.2288',
    'time': '5:19:20.92138', 'psrname': '1133+16', 'N used channels': '450',
    'date': '1/1/2017', 'period': '1.18780828035'}

    >> print(len(t_point))
    570

    >> print(len(i_point))
    570
    """

    header = read_header(filename)
    observs = np.genfromtxt(filename, skip_header=14).T
    time_point = observs[0]
    intens_point = observs[1]

    return header, time_point, intens_point


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
    delay : numpy.float64
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

    return np.round(delay, 8)


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
    TB_sec : numpy.float64
            The moment of start the recording of a given pulsar in
            the barycentric time recalculated to the barycenter
            of the solar system for a given Julian date

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


def get_time_delay_full(file_obs_1, file_obs_2, file_out_1, file_out_2, tz1=3, tz2=3):
    """
    Help on function get_time_delay_full in module gpclock:

    get_time_delay_full(file_obs_1, file_obs_2, file_out_1, file_out_2):
        return delay

    Discription
    ----------

    Parameters
    ----------
    file_obs_1 : str
        Input data. Name of file of observation in current directory
        or path to the file.

    file_obs_2 : str
        Input data. Name of file of observation in current directory
        or path to the file.

    file_out_1 : str
        Input data. Name of file of baricentric times start in current directory
        or path to the file.

    file_out_2 : str
        Input data. Name of file of baricentric times start in current directory
        or path to the file.

    tz1 : int
        Input data. Time zone for first recorder. 3 by default.

    tz2 : int
        Input data. Time zone for second recorder. 3 by default.

    Returns
    -------
    delay : numpy.float64
            Time delay between two arrays.

    Examples
    --------
    """

    header_1, _, flux_1 = read_prf(file_obs_1)
    header_2, _, flux_2 = read_prf(file_obs_2)

    day_1, month_1, year_1 = header_1['date'].split('/')
    hour_1, minute_1, second_1 = header_1['time'].split(':')
    second_1, microsecond_1 = second_1.split('.')
    time_start_1 = datetime(
        int(year_1),
        int(month_1),
        int(day_1),
        int(hour_1),
        int(minute_1),
        int(second_1),
        int(microsecond_1))
    time_start_1 -= dt.timedelta(hours=tz1)  # to UTC time
    fs_p_1, _ = str(time_start_1.to_mjd()).split('.')

    day_2, month_2, year_2 = header_2['date'].split('/')
    hour_2, minute_2, second_2 = header_2['time'].split(':')
    second_2, microsecond_2 = second_2.split('.')
    time_start_2 = datetime(
        int(year_2),
        int(month_2),
        int(day_2),
        int(hour_2),
        int(minute_2),
        int(second_2),
        int(microsecond_2))
    time_start_2 -= dt.timedelta(hours=tz2)  # to UTC time
    fs_p_2, _ = str(time_start_2.to_mjd()).split('.')

    tay = np.float64(header_1['tau'])/1000.

    ts_1 = get_TB_sec(file_out_1, fs_p_1)
    ts_2 = get_TB_sec(file_out_2, fs_p_2)

    return get_time_delay(ts_1, ts_2, flux_1, flux_2, tay)
