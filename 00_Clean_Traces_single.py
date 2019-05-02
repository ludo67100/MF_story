# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 16:29:25 2019

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

window = [0.4,0.6] #in seconds
baseline = [0.4,0.45] #in seconds
stim = [0.50,0.010] # start and duration in seconds

holdings = [-70,-50,0]


savefig = False
savedata = True
savedir = 'U:/01_ANALYSIS/MF_BINDA'

bl_start = 0

url = 'U:/RAW DATA/MF_BINDA/Sorted/19-04-2019/Cell_2'
name = '190419_2'

RAW_SIGS, AVG_SIGS = [],[] #For raw sigs as matrix and avg sigs as arrays

#-------------------------------------FUNCTIONS-------------------------------------------------------

#------------------------------------FIRST COLLECT THE DATA--------------------------------
for holding, idx in zip(holdings,range(len(holdings))):
    
    
    path = '{}/{}'.format(url,str(holding))
    file_list = sorted(os.listdir(path))
    
    print (file_list)
    
    #Basic info 
    r = neo.io.WinWcpIO('{}/{}'.format(path,file_list[0]))
    
    b=r.read_block()
    
    sampling_rate = float(b.segments[0].analogsignals[0].sampling_rate)
    time_vector = np.ravel(b.segments[0].analogsignals[0].times)
    
    
    signals = np.zeros((len(b.segments),len(time_vector))) #The matrix for raw traces
    
    I_tags = []

    
    for sweep in range(len(b.segments)):
        
        sig = np.ravel(b.segments[sweep].analogsignals[ch].magnitude)
        
        bl_begin = np.where(time_vector>=baseline[0])[0][0]
        bl_end = np.where(time_vector>=baseline[1])[0][0]
        leak = np.mean(sig[bl_begin:bl_end])
        
        signal = sig - leak

        if holding == 0 : #So inhibition, traces have to be cleaned manually 
            
            print('')
            
            print ('H=0mV, sweep#{}'.format(sweep))
            
            BL = np.nanmean(signal[bl_begin:bl_end])
                   
            plt.figure()
            plt.axvspan(stim[0],stim[0]+stim[1],color='skyblue',alpha=0.5) #For the stim
            plt.xlim(0.45,0.55)
            plt.ylim(-200,200)
            plt.plot(time_vector,np.ones(len(time_vector))*BL,label='baseline')
            plt.plot(time_vector,signal,label='signal')
            plt.legend(loc='best')
            plt.pause(0.1)
            plt.show(block=False)
        
            decide = input('Keep this sweep ? : y/n : ')
            
            if decide=='break':
                sys.exit()
            
            if decide == 'n':      
                nans = np.zeros(len(signal))
                nans[:] = np.nan
                signals[sweep,:] = nans
                I_tags.append(0)
                plt.close()
                
            else:
                signals[sweep,:] = signal
                I_tags.append(1)
                plt.close()

        signals[sweep,:] = signal
     
    #Append traces to matrices for later
    RAW_SIGS.append(signals)

print ('MANUAL SORTING: DONE.')
print ('Computing & saving...')

    
I_tags_df = pd.DataFrame(I_tags)
I_tags_df.to_excel('{}/{}_inhibition_tag_list.xlsx'.format(savedir,name))

#------------------------------Now the serious shit----------------------------

with pd.ExcelWriter('{}/{}_cleaned_avg_recordings.xlsx'.format(savedir,name)) as writer:
        
    for holding, holding_idx in zip(holdings, range(len(holdings))):
    
        AVG_SIGS = np.nanmean(RAW_SIGS[holding_idx],axis=0)
        plt.figure()
        plt.plot(time_vector,AVG_SIGS)
        plt.title('{} : average recording at H={}mV'.format(name,holding))

        #Make col names for dataframe and index 
        indexes = time_vector

        #Put the sigs in a dataframe
        Average_sigs = pd.DataFrame(np.asarray(AVG_SIGS).transpose(),
                                    index=indexes)
        
        
        #Save the data (or not)
        if savedata == True:
            Average_sigs.to_excel(writer,sheet_name='H={}mV'.format(holding),na_rep='nan')
























    