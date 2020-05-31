# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 18:08:43 2020

@author: apan
"""

from char_compare import Manager
from char_compare import data_slicer as ds
from char_compare import langmuir_models2 as lm
from Fluctuation_analysis import lp_spectrum as ls
import numpy as np
import matplotlib.pyplot as plt

def running_mean(x, N):
        N = int(N)
        #y_array = np.array([])
        ll = int(round(len(x)/N))
        y_array = np.zeros(shape = ll)#, dtype=complex)
        
        for l in range (0,ll):
            y = np.sum(x[N*l:N*(l+1)])/N
            y_array[l] =  y
        return (y_array)

def create_series(blp):
        blp.voltage *= 100
        blp.current *= -2
        vmcs = ds.VmaxClusterSlicer(blp, offset = -40, half_periods = 2)
        model = lm.WeinlichSheathExpansion
        blp.fit(model, vmcs)
        
        Te = blp.fit_result(model, vmcs, 'params', 'Te')
        ne = blp.fit_result(model, vmcs, 'params', 'ne')
        Vf = blp.fit_result(model, vmcs, 'params', 'Vf')
        time = blp.fit_result(model, vmcs, 'time')
        Vp = Vf + 3.0*(Te)
        
        return (time, ne, Te, Vf, Vp)
               

def cross(time, fseries, gseries, nave):
        l = len(fseries)
        Fw = np.zeros(shape = (l), dtype=complex)
        Gw = np.zeros(shape = (l), dtype=complex)
        Pfg = np.zeros(shape = (l), dtype=complex)
        Pff = np.zeros(shape = (l), dtype=complex)
        Pgg = np.zeros(shape = (l), dtype=complex) 
        freq = np.linspace(0, 1, round(len(time)/2))*(0.5/(time[2]-time[1]))
        freq = ls.running_mean(freq, nave)
        Fw = (np.fft.fft(fseries-np.average(fseries)))
        Gw = (np.fft.fft(gseries-np.average(gseries)))
        Pff = (np.conj(Fw)*Fw)
        Pgg = (np.conj(Gw)*Gw)
        Pfg = (np.conj(Fw)*Gw)
        cfg = (abs(Pfg)/((Pff*Pgg)**0.5))
        cross_phase = (np.angle(Pfg))
        ##for plotting
        param1 = running_mean(Fw, nave)
        param2 = running_mean(Gw, nave)
        param3 = running_mean(Pfg, nave)
        param4 = running_mean(cfg, nave)
        
        fig, axs = plt.subplots(3, 2, figsize=(12, 10))
        
        axs[0, 0].semilogy(freq, abs(param1[:len(freq)])/max(abs(param1)))
        axs[0, 0].set_title('FFT of signal 1')
        axs[0, 0].set(xlabel='frequency (Hz)', ylabel='Power density (a.u.)')
        
        axs[0, 1].semilogy(freq, abs(param2[:len(freq)])/max(abs(param2)), 'tab:orange')
        axs[0, 1].set_title('FFT of signal 2')
        axs[0, 1].set(xlabel='frequency (Hz)', ylabel='Power density (a.u.)')
        
        axs[1, 0].semilogy(freq, abs(param3[:len(freq)])/max(abs(param3)), 'tab:green')
        axs[1, 0].set(xlabel='frequency (Hz)', ylabel='Cross power density (a.u.)')
        
        axs[1, 1].semilogy(freq, param4[:len(freq)], 'tab:red')
        axs[1, 1].set(xlabel='frequency (Hz)', ylabel='Coherence')
        
        axs[2, 0].plot(freq, cross_phase[:len(freq)]*180/np.pi)
        axs[2, 0].set(xlabel ='frequency (Hz)', ylabel='Cross Phase (degrees)')
        
        return (freq, param4, param1)
    
    
def covariance(xseries, yseries):
        #xmean, ymean = np.mean(xseries), np.mean(yseries)
        xsd, ysd = xseries-np.mean(xseries), yseries - np.mean(yseries)
        cov = np.mean(xsd*ysd)
        return (cov)
    
def autocovariance(dt, xseries, yseries, nmax, n):
        cov = np.zeros(shape = nmax)
        for i in range (nmax):
            xseries1 = xseries[0:len(xseries)-n*i]
            yseries1 = yseries[i*n:]
            #print (i, len(xseries1), len(yseries1))
            cov[i] = covariance(xseries1, yseries1)
        return (cov)

def autocovariance2(Xt):
        n = int(len(Xt))
        Xmu = np.zeros(shape = n)
        for m in range (n):
            imu = 0
            for u in range (int(round((n+1)/(m+1)))):
                imu += Xt[u*m-m+1]
            Xmu[m] = imu/(m+1)
        return (Xmu)
    
def findpdf(xseries, index, no):
        dur = int(len(xseries))
        fmax = max(abs((xseries)-np.mean(xseries)))
        rnge  = fmax/(index)
        Pa = np.zeros(shape = 2*index)
        sigma = np.zeros(shape = 2*index)
        for i in range (2*index):
            xa = np.mean(xseries)-fmax+rnge*(i)
            xb = np.mean(xseries)-fmax+rnge*(i+1)
            p = 0
            for l in range (dur):
                x = xseries[l]
                if x>=xa and x<xb:
                    p +=1
            sigma[i] = xa-np.mean(xseries)
            Pa[i] = p
            #print (i, l)
        prob = Pa/sum(Pa)
        prob = running_mean(prob, no)
        sigma = running_mean(sigma, no)
        skewness = np.mean((prob-np.mean(prob))**3)/(np.mean((prob-np.mean(prob))**2))**1.5
        kurtosis = np.mean((prob-np.mean(prob))**4)/(np.mean((prob-np.mean(prob))**2))**2
        print('Skewness = {}, Kurtosis = {}'.format(skewness, kurtosis))
        plt.plot(sigma[10:-10], prob[10:-10], 'b')
        plt.xlabel(r'$n-\~n$')
        plt.ylabel('PDF')
        plt.show()
                
"""sdef cross_spec(tseries, fseries, gseries, bin_time, npoints):
        farray = fseries#ls.bins(tseries, fseries, bin_time)
        garray = gseries#ls.bins(tseries, gseries, bin_time)
        time = tseries#ls.bins(tseries, tseries, bin_time)

        #stime = time_array[0]
        l = len(farray)
        m = len(farray[0])
        Fw = np.zeros(shape = (l, m), dtype=complex)
        Gw = np.zeros(shape = (l, m), dtype=complex)
        Pfg = np.zeros(shape = (l, m), dtype=complex)
        Pff = np.zeros(shape = (l, m), dtype=complex)
        Pgg = np.zeros(shape = (l, m), dtype=complex)
        
        freq = np.linspace(0, 1, round(len(time)/2))*(0.5/(time[2]-time[1]))    
        nave = len(freq)/npoints
        for i  in range (0, l):
            Fw[i] = np.fft.fft(farray[i])
            Gw[i] = np.fft.fft(garray[i])
            Pff[i] = np.conj(Fw[i])*Fw[i]
            Pgg[i] = np.conj(Gw[i])*Gw[i]
            Pfg[i] = np.conj(Fw[i])*Gw[i]
        freq = ls.running_mean(freq, nave)    
        nFw = np.zeros(shape = m, dtype=complex)
        nGw = np.zeros(shape = m, dtype=complex)
        nPfg = np.zeros(shape = m, dtype=complex)
        Fw_ave = np.zeros(shape = m, dtype=complex)
        Gw_ave = np.zeros(shape = m, dtype=complex)
        Pfg_ave = np.zeros(shape = m, dtype=complex)
        
        nFw = ls.array_avg(Fw)
        Fw_ave = ls.running_mean(nFw, nave)
        
        nGw = ls.array_avg(Gw)
        Gw_ave = ls.running_mean(nGw, nave)
        
        nPfg = ls.array_avg(Pfg)
        cross_phase = np.angle(nPfg)
        Pfg_ave = ls.running_mean(nPfg, nave)
        cross_phase = ls.running_mean(cross_phase, nave)
        cp = np.angle(Pfg_ave)
        
        nPff = ls.array_avg(Pff)
        nPgg = ls.array_avg(Pgg)
        Cfg = abs(nPfg)/((nPff*nPgg)**0.5)
        Cfg_ave = ls.running_mean(Cfg, nave)
        print (len(freq), len(Fw_ave))
        
        fig, axs = plt.subplots(3, 2, figsize=(20, 10))
        
        axs[0, 0].semilogy(freq, abs(Fw_ave)[:len(freq)])
        axs[0, 0].set_title('FFT of signal 1')
        axs[0, 0].set(xlabel='frequency (Hz)', ylabel='Power density (a.u.)')
        
        axs[0, 1].semilogy(freq, abs(Gw_ave)[:len(freq)], 'tab:orange')
        axs[0, 1].set_title('FFT of signal 2')
        axs[0, 1].set(xlabel='frequency (Hz)', ylabel='Power density (a.u.)')
        
        axs[1, 0].semilogy(freq, abs(Pfg_ave[:len(freq)]), 'tab:green')
        axs[1, 0].set(xlabel='frequency (Hz)', ylabel='Cross power density (a.u.)')
        
        axs[1, 1].semilogy(freq, Cfg_ave[:len(freq)], 'tab:red')
        axs[1, 1].set(xlabel='frequency (Hz)', ylabel='Coherence')
        
        axs[2, 0].plot(freq, cross_phase[:len(freq)])
        axs[2, 0].set(xlabel ='frequency (Hz)', ylabel='Cross Phase (radians)')
        """