# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 17:25:01 2020

@author: apan
"""

#BLP data analysis based on the document by Peter Kallmeyer 1-GDI00-T0106

import numpy as np
import matplotlib.pyplot as plt
#from scipy import signal
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

def running_mean(x, N):
        N = int(N)
        y_array = np.array([])
        for l in range (0,int(len(x)/N)):
            y = np.sum(x[N*l:N*(l+1)])/N
            y_array = np.append(y_array, y)
        return (y_array)

def array_avg(array):
        n = len(array)
        l = len(array[1])
        a = np.zeros(l)
        
        for i in range (0, l):
            sum = 0
            #print (i)
            for k in range (0, n):
                #print (k, i)
                #print(array[k][i])
                sum = sum + array[k][i]
                #print(sum, k, i)
            #print(sum, k, i)
            a[i] = sum/n
        return (a)

def bins(tseries, fseries, bin_time):
        ns = len(tseries) #length of series  
        dt = (tseries[2]-tseries[1])
        total_time = ns*dt # total time for which data is available
        #bin_time = 50e-3 # time interval chosen for 1 data bin
        bin_number = (total_time//bin_time) # number of bins available for total time
        #bin_points = int(bin_time/(tseries[2]-tseries[1])) # number of points in one bin
        index = int(round(bin_number*bin_time/dt))
        fseries = fseries[:index]
        print (index)
        new_matrix = np.split(fseries, bin_number) 
        return (new_matrix)
    
def blp_qty(time, Udummy, Udiff, tau = 0.3e-6):
        Rs = 50 #ohm
        dt = time[2]-time[1]
        lag = int(round(tau/dt))
        Us = Udummy - Udiff
        iprobe = (Udiff[lag:])/Rs
        Tcyc = 2e-5
        Ncyc = int(5)
        Npoints = int(round(Tcyc/dt)) #number of points in one sweep cycle (50 kHz)
        #Udummy = ((Udummy+10)*10) #let this be the scaling
        Uoffset = (max(Udummy[:Npoints])+min(Udummy[:Npoints]))/2
        Usource = ((Udummy[lag:]-Uoffset)/np.cos(2*np.pi*50000*tau))+Uoffset
        Udummy = Udummy[:Ncyc*Npoints]
        #Udummy = Udummy[lag:]
        iprobe = iprobe[:Ncyc*Npoints]
        Usource = Usource[:Ncyc*Npoints]
        
        return (Udummy, iprobe, Usource)

def data_prep(time, voltage, current, time_bin):
        time_matrix = bins(time, time, time_bin)
        voltage_matrix = bins(time, voltage, time_bin)
        current_matrix = bins(time, current, time_bin)
        
        time = time_matrix[0]
        voltage_unscaled = array_avg(voltage_matrix)
        current_unscaled = array_avg(current_matrix)
        
        voltage = data_scale(voltage_unscaled)
        voltage = voltage[:-3]
        current = data_scale(current_unscaled)
        current = -current[3:]
        return (time, voltage, current)

def data_scale(data):
        #x = min(data)
        #y = -80-100*x
        data = (data*100)
        return(data)
    
def lsfit(time, series):#, p1, p2, p3):
        
        t = running_mean(time, 10)
        data = running_mean(series, 10)
        
        def test_fit(x, a, b, c):
            return a+ b*np.sin((2*np.pi*50000*x)+c)
        
        guess_mean = np.mean(data)
        guess_phase = 1
        guess_amp = (max(data)-guess_mean)/2
        #est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, [guess_amp, guess_freq, guess_phase, guess_mean])[0]
        params, params_covariance = curve_fit(test_fit, t, data, p0=[guess_mean, guess_amp, guess_phase])
      
        data_fit = params[0] +params[1]*np.sin((2*np.pi*50000*time)+params[2])
        plt.plot(t, data, '.', label ='averaged')
        plt.plot(time,series, label='original data')
        plt.plot(time, data_fit, label='after fitting')
        plt.legend()
        plt.show()