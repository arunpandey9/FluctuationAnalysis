# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:36:57 2020

@author: apan
"""
import numpy as np
import matplotlib.pyplot as plt 
from scipy import signal


def running_mean(x, N):
        N = int(N)
        #y_array = np.array([])
        y_array = np.zeros(shape = int(round(len(x)/N)), dtype=complex)
        
        for l in range (0,int(len(x)/N)):
            y = np.sum(x[N*l:N*(l+1)])/N
            y_array = np.append(y_array, y)
        return (y_array)
    
        """ cumsum = np.cumsum(np.insert(x, 0, 0)) 
            return (cumsum[N:] - cumsum[:-N]) / float(N)
        """
        
def spectran(time, fseries, gseries, npoints, nfreq):
        freq = np.zeros(npoints)
        #f_limit = 250000/n_ave
        fseries = fseries-np.mean(fseries)
        gseries = gseries-np.mean(gseries)
        #print ('frequency limit is {}'. format(f_limit))
        #time = running_mean(time, n_ave)
        fsampling = 1/(time[2]-time[1])
        #total_time = len(fseries)/fsampling
        #fseries = running_mean(fseries, n_ave)
        #gseries = running_mean(gseries, n_ave)
        freq_1, Pxy = signal.csd(fseries, gseries, fs = fsampling, nperseg = npoints)
        freq, cxy = signal.coherence(fseries, gseries, fs = fsampling, nperseg = npoints)
        Fw = np.fft.fft(fseries)
        Gw = np.fft.fft(fseries)
        freq_fft = np.linspace(0, 1, int(len(time)/2))*(0.5/(time[2]-time[1]))
        
        fig, axs = plt.subplots(2, 2)
        n_ind = (int(len(Fw)/npoints))
        freq_fft = running_mean(freq_fft, n_ind*nfreq)
        freq = running_mean(freq, nfreq)
        Fw = running_mean(Fw, n_ind*nfreq)
        Gw = running_mean(Gw, n_ind*nfreq)
        Pxy = running_mean(Pxy, nfreq)
        cxy = running_mean(cxy, nfreq)
        print(len(Fw), len(Pxy))
        axs[0, 0].semilogy(freq_fft[:int(len(Fw)/2)], abs(Fw)[:int(len(Fw)/2)])
        axs[0, 0].set_title('FFT of signal 1')
        #axs[0, 0].set_ylabel('Power density (a.u.)')
        #axs[0, 0].set_xlabel('frequency (Hz)')
        axs[0, 0].set(xlabel='frequency (Hz)', ylabel='Power density (a.u.)')
        
        axs[0, 1].semilogy(freq_fft[:int(len(Fw)/2)], abs(Gw)[:int(len(Fw)/2)], 'tab:orange')
        axs[0, 1].set_title('FFT of signal 2')
        #axs[0, 1].ylabel('Power density (a.u.)')
        #axs[0, 1].xlabel('frequency (Hz)')
        axs[0, 1].set(xlabel='frequency (Hz)', ylabel='Power density (a.u.)')
        
        axs[1, 0].semilogy(freq, abs(Pxy), 'tab:green')
        #axs[1, 0].set_title('Cross spectrum density')
        #axs[1, 0].ylabel('Cross spectrum density (a.u.)')
        #axs[1, 0].xlabel('frequency (Hz)')
        axs[1, 0].set(xlabel='frequency (Hz)', ylabel='Cross spectrum density (a.u.)')
        
        axs[1, 1].semilogy(freq, cxy, 'tab:red')
        #axs[1, 1].set_title('Coherence]')
        #axs[1, 1].ylabel('Coherence')
        #axs[1, 1].xlabel('frequency (Hz)')
        axs[1, 1].set(xlabel='frequency (Hz)', ylabel='Coherence')
        """for ax in axs.flat:
            ax.set(xlabel='x-label', ylabel='y-label')

        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()
        
        plt.figure()
        plt.semilogy(freq, abs(Pxy))
        plt.ylabel('Cross spectrum density (a.u.)')
        plt.xlabel('frequency (Hz)')
        plt.show()
        
        plt.figure()
        plt.ylim(0.0001, 2)
        plt.semilogy(freq, cxy)
        plt.ylabel('Coherence')
        plt.xlabel('frequency (Hz)')
        plt.show()
        """
        return freq, Pxy, cxy

def datasort(tseries, yseries, f_limit):
        N = int(250000/f_limit) # no. of points on which moving average is performed 
        print("averaging over {} points".format(N)) 
        ns = len(tseries)   
        nt = ns*(tseries[2]-tseries[1]) # total time for which data is available
        dt1 = 50e-3 # time interval chosen for 1 data bin
        dn1 = int(nt/dt1) # number of bins available for total time
        npoints = int(dt1/(tseries[2]-tseries[1])) # number of points in one bin
        for i in range (0, dn1):
            time = tseries[npoints*i:npoints*(i+1)]
            current = yseries[npoints*i:npoints*(i+1)]
            current = current -np.mean(current)
            time_avg = running_mean(time, N)
            #freq = np.linspace(0,1, len(time))*(1/(time[2]-time[1]))
            current_avg = running_mean(current, N)
            current_fft_complex = (np.fft.fft(current_avg))
            current_fft = abs(current_fft_complex)
            Na = len(time_avg)
            #dt = (time_avg[2]-time_avg[1])
            print(i+1)
            freq = np.linspace(0, 1, int(Na/2))*(0.5/(time_avg[2]-time_avg[1]))   
            plt.figure()
            plt.plot(freq, current_fft[:int(Na/2)])
            # Show/save figure as desired.
            plt.show()
            
            #print (time[1:10], len(time))
            #return (current_fft_complex)

def fastft(tseries, fseries):
        fseries = fseries - np.mean(fseries)
        Fw = np.fft.fft(fseries)
        freq = np.linspace(0, 1, round(len(tseries)/2))*(0.5/(tseries[2]-tseries[1]))
        return (freq, Fw)
    

def bins(tseries, fseries, bin_time):
        ns = len(tseries) #length of series  
        total_time = ns*(tseries[2]-tseries[1]) # total time for which data is available
        #bin_time = 50e-3 # time interval chosen for 1 data bin
        bin_number = round(total_time/bin_time) # number of bins available for total time
        #bin_points = int(bin_time/(tseries[2]-tseries[1])) # number of points in one bin
        new_matrix = np.split(fseries, bin_number) 
        return (new_matrix)
    
def cross_spec(tseries, fseries, gseries, bin_time, npoints):
        farray = bins(tseries, fseries, bin_time)
        garray = bins(tseries, gseries, bin_time)
        time_array = bins(tseries, tseries, bin_time)
        time = time_array[0]
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
        freq = running_mean(freq, nave)    
        nFw = np.zeros(shape = m, dtype=complex)
        nGw = np.zeros(shape = m, dtype=complex)
        nPfg = np.zeros(shape = m, dtype=complex)
        Fw_ave = np.zeros(shape = m, dtype=complex)
        Gw_ave = np.zeros(shape = m, dtype=complex)
        Pfg_ave = np.zeros(shape = m, dtype=complex)
        
        nFw = array_avg(Fw)
        Fw_ave = running_mean(nFw, nave)
        
        nGw = array_avg(Gw)
        Gw_ave = running_mean(nGw, nave)
        
        nPfg = array_avg(Pfg)
        cross_phase = np.angle(nPfg)
        Pfg_ave = running_mean(nPfg, nave)
        cross_phase = running_mean(cross_phase, nave)
        cp = np.angle(Pfg_ave)
        
        nPff = array_avg(Pff)
        nPgg = array_avg(Pgg)
        Cfg = abs(nPfg)/((nPff*nPgg)**0.5)
        Cfg_ave = running_mean(Cfg, nave)
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
        
        #axs[2, 1].plot(freq, cp[:len(freq)])
        #axs[2, 1].set(xlabel ='frequency (Hz)', ylabel='Cross Phase (radians)')
        plt.savefig('figure.pdf')
        #return (freq, Pfg_ave, nPfg, Cfg_ave)

def array_avg(array):
        n = len(array)
        l = len(array[1])
        a = np.zeros(l, dtype = complex)
        
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

"""def cross_spec(tseries, fseries, gseries, f_limit):
        N = int(250000/f_limit) # no. of points on which moving average is performed 
        print("averaging over {} points".format(N)) 
        ns = len(tseries)   
        total_time = ns*(tseries[2]-tseries[1]) # total time for which data is available
        time_bin = 50e-3 # time interval chosen for 1 data bin
        num_bins = int(total_time/time_bin) # number of bins available for total time
        npoints = int(time_bin/(tseries[2]-tseries[1])) # number of points in one bin
        for i in range (0, num_bins):
            time = tseries[npoints*i:npoints*(i+1)]
            current_f = fseries[npoints*i:npoints*(i+1)]
            current_f = current_f -np.mean(current_f)
            current_g = gseries[npoints*i:npoints*(i+1)]
            current_g = current_g -np.mean(current_g)
            time_avg = running_mean(time, N)
            Na = len(time_avg)
            #dt = (time_avg[2]-time_avg[1])
            print(i+1)
            freq = np.linspace(0, 1, int(Na/2))*(0.5/(time_avg[2]-time_avg[1]))
            #freq = np.linspace(0,1, len(time))*(1/(time[2]-time[1]))
            current_avg_f = running_mean(current_f, N)
            current_avg_g = running_mean(current_g, N)
            current_fft_complex_f = (np.fft.fft(current_avg_f))
            current_fft_complex_g = (np.fft.fft(current_avg_g))
            P_fg = (np.conj(current_fft_complex_f))*(current_fft_complex_g)
            P_abs = abs(P_fg)#P_abs = abs(current_fft_complex_1)*abs(current_fft_complex_2)
            #C_fg = np.real(P_fg)
            #Q_fg = np.imag(P_fg)
            P_ff = (np.conj(current_fft_complex_f))*(current_fft_complex_f)
            P_gg = (np.conj(current_fft_complex_g))*(current_fft_complex_g)
            gamma_fg = P_abs/(P_ff*P_gg)**0.5
            cross_spect = abs(gamma_fg)
            #print(P_abs[i], P_ff[i], P_gg[i], gamma_fg[i])
            plt.figure()
            plt.semilogy(freq, P_abs[:int(Na/2)])
            plt.ylabel('Power (a.u.)')
            plt.xlabel('frequency (Hz)')
            plt.show()

def writeinto(a, b):
        f = open('file.txt', 'w')
        for i in range(len(a)):
            f.write("%f %f\n" % (a[i], b[i]))
        f.close()"""