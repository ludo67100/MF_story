# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 17:06:33 2019

Computes figures for group data in MF experiments for Binda et al 

@author: ludovic.spaeth
"""

import numpy as np
import pandas as pd 
from matplotlib import pyplot as plt 

#-----------------------------------FUNCTIONS--------------------------------------------
def LinReg(x,y,conf=0.95,plot=True):
    '''
    Computes linear regression on 2 arrays 
    
    x,y (1D arrays) = the x & y data 
    
    conf (float) = confidence treshold (alpha = 1-conf)
    
    plot (bool) : if True, will plot the result of the linear reg 
    
    Returns : 
        
        px (array) : x-axis for plot, based on min/max values of x
        
        nom (array) : y-values of linear model (nomial values of y)
        
        lpb, upb (arrays) : lower and upper prediction bands 
        
        r2 = squared correlation coefficient 
    
    '''
    import numpy as np 
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    from scipy import stats
    import uncertainties as unc
    import matplotlib
    from sklearn import metrics
    matplotlib.rcParams['pdf.fonttype'] = 42
    
    #Pip install uncertainties if needed 
    try :
        import uncertainties.unumpy as unp
        
    except : 
        import pip
        pip.main(['install','uncertainties'])
        import uncertainties.unumpy as unp 
    
    
    n = len(y)
    
    def f(x, a, b):
        return a * x + b
    
    popt, pcov = curve_fit(f, x, y)
    
    # retrieve parameter values
    a = popt[0]
    b = popt[1]
    print('Optimal Values')
    print('a: ' + str(a))
    print('b: ' + str(b))
    
    # compute r^2
    r2 = 1.0-(sum((y-f(x,a,b))**2)/((n-1.0)*np.var(y,ddof=1)))
    print('R^2: ' + str(r2))
    
    # calculate parameter confidence interval
    a,b = unc.correlated_values(popt, pcov)
    print('Uncertainty')
    print('a: ' + str(a))
    print('b: ' + str(b))
    
    # plot data
    if plot == True:
        plt.figure()
        plt.scatter(x, y, s=20, label='Data')
    
    # calculate regression confidence interval
    px = np.linspace(np.min(x),np.max(x),n)
    py = a*px+b
    nom = unp.nominal_values(py)
    std = unp.std_devs(py)
    
    def predband(x, xd, yd, p, func, conf=conf):
        '''
        x = requested points
        xd = x data
        yd = y data
        p = parameters
        func = function name
        '''
        alpha = 1.0 - conf    # significance
        N = xd.size          # data sample size
        var_n = len(p)  # number of parameters
        # Quantile of Student's t distribution for p=(1-alpha/2)
        q = stats.t.ppf(1.0 - alpha / 2.0, N - var_n)
        # Stdev of an individual measurement
        se = np.sqrt(1. / (N - var_n) * \
                     np.sum((yd - func(xd, *p)) ** 2))
        # Auxiliary definitions
        sx = (x - xd.mean()) ** 2
        sxd = np.sum((xd - xd.mean()) ** 2)
        # Predicted values (best-fit model)
        yp = func(x, *p)
        # Prediction band
        dy = q * se * np.sqrt(1.0+ (1.0/N) + (sx/sxd))
        # Upper & lower prediction bands.
        lpb, upb = yp - dy, yp + dy
        return lpb, upb
    
    lpb, upb = predband(px, x, y, popt, f, conf=conf)
        
    if plot == True:
        # plot the regression
        plt.plot(px, nom, c='orange', label='y=a x + b',linewidth=2)
        
        # uncertainty lines (95% confidence)
        plt.plot(px, nom - 1.96 * std, c='0.5',linestyle='--',\
                 label='95% Confidence Region')
        plt.plot(px, nom + 1.96 * std, c='0.5',linestyle='--')
        # prediction band (95% confidence)
        plt.plot(px, lpb, color='0.5',label='95% Prediction Band',linestyle=':')
        plt.plot(px, upb, color='0.5',linestyle=':')
        plt.ylabel('Y')
        plt.xlabel('X')
        plt.legend(loc='best')
        plt.title('Linear Reg : R$^2$={}'.format(round(r2,2)))
        plt.show()
        
    return px,nom,lpb,upb,r2,std


#-------------------------------------CODE---------------------------------------

url = 'U:/01_ANALYSIS/MF_BINDA/00_Half_Manual_latency_mesures_MF_BINDA_Ludo_DATASET.xlsx'

cl_datasheet = pd.read_excel(url, sheet_name='hard clean')


delta_lat = cl_datasheet[['D_lat (ms)']].values.ravel()
EPSQ = cl_datasheet[['EPSQ(pC)']].values.ravel()
IPSQ = cl_datasheet[['IPSQ(pC)']].values.ravel()

real_ffi, outsiders = [],[]

in_EPSQ, in_IPSQ = [],[]
out_EPSQ, out_IPSQ = [],[]

for i in range(len(delta_lat)):
    if 0. <= delta_lat[i] <= 5.:
        real_ffi.append(delta_lat[i])
        in_EPSQ.append(np.abs(EPSQ[i]))
        in_IPSQ.append(IPSQ[i])
        
    else:
        outsiders.append(delta_lat[i])
        out_EPSQ.append(np.abs(EPSQ[i]))
        out_IPSQ.append(IPSQ[i])



fig, ax = plt.subplots(1,3,figsize=(13,5))
fig.suptitle('MF Stimulation - cleaned recordings')

#Delta lat histogramm
bins = np.arange(-15,15,0.5)
ax[0].hist(real_ffi, bins=bins, label='group 1', color='green')
ax[0].hist(outsiders, bins=bins, label='group 2', color='lightcoral')
ax[0].legend(loc='best')
ax[0].set_xlabel(r'$\Delta$ latency (ms)'), ax[0].set_ylabel('Count')

#Group 1 scatter
ax[1].scatter(in_EPSQ, in_IPSQ,color='green',label='group 1')
px,nom,lpb,upb,r2,std = LinReg(np.asarray(in_EPSQ),np.asarray(in_IPSQ),conf=0.95,plot=False)
ax[1].plot(px,nom,color='green',label='y=ax+b',linewidth=2)
ax[1].plot(px, nom - 1.96 * std, c='0.5',linestyle='--',\
         label='95% Confidence Region')
ax[1].plot(px, nom + 1.96 * std, c='0.5',linestyle='--')
ax[1].plot(px, lpb, color='0.5',label='95% Prediction Band',linestyle=':')
ax[1].plot(px, upb, color='0.5',linestyle=':')
ax[1].set_xlabel('EPSQ (pC)'); ax[1].set_ylabel('IPSQ (pC)')
ax[1].set_title('R={} / $R^2$={}'.format(round(np.sqrt(r2),2),round(r2,2)))

#Group 2 scatter
ax[2].scatter(out_EPSQ, out_IPSQ,color='lightcoral',label='group 2')
px,nom,lpb,upb,r2,std = LinReg(np.asarray(out_EPSQ),np.asarray(out_IPSQ),conf=0.95,plot=False)
ax[2].plot(px,nom,color='lightcoral',label='y=ax+b',linewidth=2)
ax[2].plot(px, nom - 1.96 * std, c='0.5',linestyle='--',\
         label='95% Confidence Region')
ax[2].plot(px, nom + 1.96 * std, c='0.5',linestyle='--')
ax[2].plot(px, lpb, color='0.5',label='95% Prediction Band',linestyle=':')
ax[2].plot(px, upb, color='0.5',linestyle=':')
ax[2].set_xlabel('EPSQ (pC)'); ax[2].set_ylabel('IPSQ (pC)')
ax[2].set_title('R={} / $R^2$={}'.format(round(np.sqrt(r2),2),round(r2,2)))