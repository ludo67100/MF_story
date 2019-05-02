# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 16:31:05 2019

@author: ludovic.spaeth
"""

import numpy as np
import neo 
import os 
import pandas as pd 
from matplotlib import pyplot as plt 
import sys

grid = np.array([[1,2,3,4],[5,6,7,8]])

ch = 0

window = [0.5,0.6] #in seconds
baseline = [0.498,0.5] #in seconds
stim = [0.50,0.010] # start and duration in seconds

holdings = [-70,-50,0]

savefig = False
savedata = True
savedir = 'U:/01_ANALYSIS/MF_BINDA/latencies'

bl_start = 0

url = 'U:/01_ANALYSIS/MF_BINDA'
cell = '180419_3b'
#----------------------------------------FUNCTIONS---------------------------------------
def latency(signal,time_vector,bl_start=0.45,bl_stop=0.49, win_begin=0.5,win_end=0.7,polarity='IPSC',
            trigger=0.5,rise_on=0.2,rise_off=0.8,extent=30,plot=True,title=None):

    '''
    Function to compute latency of PSCs based on linear fit betwwen rise_on and rise_off % of PSCs slope
    AND synaptic charge
    
    --------INPUTS--------------------------------------------------------------
    signal (1D array) =  the signal
    time_vector (1D array) = the time vector
    bl_start/bl_stop (float) = time window for baseline
    win_begin/win_end (float) = time window for stim 
    polarity (str) = 'EPSC' or 'IPSC'
    trigger (float) = time (in seconds) of stimulus start 
    rise_on (float) = % of slope to start fit 
    rise_off (float) = % of slope to stop fit
    extent (int) = n border to extent fit, then it can cross the baseline
    plot (bool) = if True, will give a plot after computation 
    
    -------OUTPUTS--------------------------------------------------------------
    latency,charge (floats) = computed latency and synaptic charge 
    
    
    '''
    import numpy as np
    from scipy.integrate import trapz 
    print ('-----------Latency computation--------------')
    
    
    win_begin = np.where(time_vector>=win_begin)[0][0]
    win_end = np.where(time_vector>=win_end)[0][0] 
    bl_begin = np.where(time_vector>=bl_start)[0][0]
    bl_end = np.where(time_vector>=bl_stop)[0][0]
    
    trunc_time_vector = time_vector[win_begin:win_end]
    
    baseline = np.mean(signal[bl_begin:bl_end],axis=0)
    print ('Average baseline : ', baseline)
    avg_bl = np.ones(len(trunc_time_vector))*baseline
    
    trunc_signal = signal[win_begin:win_end]
    charge = trapz(trunc_signal, dx = trunc_time_vector[2]-trunc_time_vector[1])
    print ('Synaptic charge = {} pC'.format(charge))

    if polarity=='EPSC':
        min_hit = np.min(trunc_signal,axis=0) #The min value = max value of the EPSC
        print ('EPSC peak : ', min_hit)
    else:
        min_hit = np.max(trunc_signal,axis=0) #The min value = max value of the EPSC
        print ('IPSC peak : ', min_hit)
        
    min_hit_index = trunc_time_vector[np.where(trunc_signal==min_hit)[0][0]] #The index of the max 
    print ('peak index : ', min_hit_index)
    
    per_80 = min_hit*rise_off
    if polarity=='EPSC':    
        per_80_index = np.where(trunc_signal<=per_80)[0][0]
    else:
        per_80_index = np.where(trunc_signal>=per_80)[0][0]
    per_80_index_time = trunc_time_vector[np.where(trunc_signal<=per_80)[0][0]] #The index of the 80% slop
    print ('80% slope coordinates : ',per_80_index_time,per_80_index)
    
    per_20 = min_hit*rise_on
    if polarity=='EPSC':    
        per_20_index = np.where(trunc_signal<=per_20)[0][0]
    else:
        per_20_index = np.where(trunc_signal>=per_20)[0][0]

    per_20_index_time = trunc_time_vector[np.where(trunc_signal<=per_20)[0][0]] #The index of the 80% slope
    print ('20% slope coordinates : ',per_20_index_time,per_20_index)
    
    slope20_80 = trunc_signal[per_20_index:per_80_index]
    time_20_80 = trunc_time_vector[per_20_index:per_80_index]
    
    from scipy.optimize import curve_fit
    
    x = time_20_80

    y = slope20_80

    def f(x, a, b):
        return a * x + b
    
    popt, pcov = curve_fit(f, x, y)
    
    # retrieve parameter values
    a = popt[0]
    b = popt[1]
    
    print ('Linear fit result : y = {}x + {}'.format(round(a,3),round(b,3)))

    ext_x = trunc_time_vector[per_20_index-extent:per_80_index+extent]
    latency = (baseline-b)/a
    print ('latency = {} s'.format(round(latency,5)))
 
    
    if plot==True:
        from matplotlib import pyplot as plt 
        plt.figure()
        plt.title(title)
        plt.plot(trunc_time_vector,trunc_signal,label='signal')
        plt.plot(trunc_time_vector,avg_bl,label='baseline')
        plt.scatter(min_hit_index,min_hit,color='r',label='slope max')
        plt.scatter(per_80_index_time,per_80,color='orange',label='slope 80%')
        plt.scatter(per_20_index_time,per_20,color='green',label='slope 20%')
        plt.plot(time_20_80,slope20_80,color='0.2',label='slope 20-80%')
        plt.plot(ext_x,f(ext_x,a,b),label='linear fit')
        plt.scatter(latency, baseline, label='Latency')
        plt.legend(loc='best')
        
    return  latency,charge


#-------------------------------------------------------------------------------

E_dataset = pd.read_excel('{}/{}_cleaned_avg_recordings.xlsx'.format(url,cell),sheet_name='H=-70mV')
I_dataset = pd.read_excel('{}/{}_cleaned_avg_recordings.xlsx'.format(url,cell),sheet_name='H=0mV')

time_vector = E_dataset.index.values

win_begin = np.where(time_vector>=window[0])[0][0]
win_compute_begin = np.where(time_vector>=0.5)[0][0]
win_end = np.where(time_vector>=window[1])[0][0] 

bl_begin = np.where(time_vector>=baseline[0])[0][0]
bl_end = np.where(time_vector>=baseline[1])[0][0]

E_LAT, I_LAT, E_CHARGE, I_CHARGE = [],[],[],[]

#Get signal
for site in range(E_dataset.shape[1]):
    plt.figure()
    plt.title('{} site#{}: I vs E'.format(cell,site))
    
    EPSC = E_dataset.iloc[:,site].values
    time = time_vector
    
    IPSC = I_dataset.iloc[:,site].values
    time = time_vector

    
    plt.plot(time[win_begin:win_end],EPSC[win_begin:win_end],label='EPSC')
    plt.plot(time[win_begin:win_end],IPSC[win_begin:win_end],label='IPSC')

    plt.axvspan(stim[0],stim[0]+stim[1],alpha=0.1,color='blue')
                                 
    E_lat, EPSQ = latency(EPSC,time,bl_start=0.49,bl_stop=0.5, win_begin=0.501,win_end=0.7,polarity='EPSC',
                          trigger=0.5,rise_on=0.1,rise_off=0.7,extent=30,plot=True,title='{} site#{}: E latency'.format(cell,site))

    I_lat, IPSQ = latency(IPSC,time,bl_start=0.49,bl_stop=0.5, win_begin=0.501,win_end=0.7,polarity='IPSC',
                          trigger=0.5,rise_on=0.1,rise_off=0.7,extent=30,plot=True,title='{} site#{}: I latency'.format(cell,site))
                                
    print('--------------------------')
    print ('site#',site,'E_lat=',E_lat)
    print ('site#',site,'I_lat=',I_lat)
    
    E_LAT.append(E_lat)
    I_LAT.append(I_lat)
    E_CHARGE.append(EPSQ)
    I_CHARGE.append(IPSQ)
    
    
    plt.legend(loc='best')
    
    
Delta_lat = np.asarray(I_LAT)-np.asarray(E_LAT)
plt.figure()
plt.hist(Delta_lat, bins=10)