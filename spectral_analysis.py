# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 12:58:22 2020

@author: apan
"""

import numpy as np
import matplotlib.pyplot as plt 

def running_mean(x, N):
        y_array = np.array([])
        for l in range (0,int(len(x)/N)):
            y = np.sum(x[N*l:N*(l+1)])/N
            y_array = np.append(y_array, y)
        return (y_array)

def fastft(tseries, yseries):
        yseries = yseries - np.mean(yseries)
        sp = (np.fft.fft(yseries))
        freq = np.linspace(0,1, len(tseries))*(1/(tseries[2]-tseries[1]))
        return (freq, sp)

def time_bin(tseries, fseries, tbin):
        total_time = len(fseries)*(tseries[2]-tseries[1])
        number_bins = int(total_time/tbin)
        npoints =  int(tbin/(tseries[2]-tseries[1])) # number of points in one bin time tbin
        current = np.array([])
        time = np.array([])
        for i in range (0, number_bins):
            time[i] = tseries[npoints*i:npoints*(i+1)]
            current[i] = fseries[npoints*i:npoints*(i+1)]
            #time[i] = np.append(time, t)
            #current[i] = np.append(current, f)
        return (time, current) 

def writeinto(a, b):
        f = open('file.txt', 'w')
        for i in range(len(a)):
            f.write("%f %f\n" % (a[i], b[i]))
        f.close()