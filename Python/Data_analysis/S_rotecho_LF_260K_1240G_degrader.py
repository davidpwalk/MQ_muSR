# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 16:29:27 2020

@author: doll_a
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy
import mufuns
from scipy import interpolate

from scipy import signal
from scipy.fftpack import fft, fftshift, fftfreq
from scipy import optimize
import iminuit
from iminuit import cost


exp_c = []

# exp_c.append(np.arange(3867,3869)) # alpha runs

# on/off runs at 1240 G with 50+50 pulse at -13dBm input
exp_c.append(np.arange(3879,3884))


run_year = 2025
# alpha taken from run 3756
run_alpha = [1.05 , 1.10 , 1.10]


for ii_exp, runs in enumerate(exp_c):

    B_vec = np.zeros_like(runs, dtype=float)

    amps = np.zeros((len(runs),3))
    stds = np.zeros((len(runs),3))
    amps_x = np.zeros((len(runs),3))

    # not used, but append_new only works with lists...
    pulspars = np.zeros((len(runs),5))
    

    for ii_run, run in enumerate(runs):

        # get the data, but actually more for the purpose of the header
        asy, asy_x, header, asy_err = mufuns.get_asy(run,run_year)

        B_vec[ii_run] = float(header['field'][:-4])

        # try to fetch info about pulse delays
        title = header['title']
        
        # try to get the pulse parameters
        try:
            pstr = title.split(',')[3]
            pdrs = pstr.split('+')
            
            pulspars[ii_run,0] = float(pdrs[0])
            pulspars[ii_run,1] = float(pdrs[1])
            
            # put the rest: 
                
            # 3rd one is direclty the sum to ease plotting
            # 4th one is indicator if it is a rotecho pulse or normal one
            
            pulspars[ii_run,2] = pulspars[ii_run,0] + pulspars[ii_run,1] 
            
            pulspars[ii_run,3] = pulspars[ii_run,1] > 0.0
            
            
            # try also the amplitude
            pulspars[ii_run,4] = -100.0 # default value: uwave off (just set to -100 dBm...)
            
            
            pstr = title.split(',')[5]
            pulspars[ii_run,4] = float(pstr.replace('dBm','').replace('m', '-'))
            
            
        except:
            pass

        asy_id = 0
        dec = 144
        
        # get the integral of the entire curve
        # integration window
        strttime = 700
        intglen = 1000
        seq_offs = 0
        
        
        firstpls_lims = np.array([strttime, strttime+intglen]) + seq_offs
        firstpls = np.where((asy_x > firstpls_lims[0]) & (asy_x < firstpls_lims[1]))[0]
        first, first_x, header, first_std = mufuns.get_asy(run, start_t = firstpls_lims[0], end_t = firstpls_lims[1], single_bin=True, year = run_year, alpha = run_alpha)
        
        strttime = 10
        intglen = 270
        seq_offs = 0
        
        # get the pre-pulse part
        prepls_lims = np.array([strttime, strttime+intglen]) + seq_offs
        prepls = np.where((asy_x > prepls_lims[0]) & (asy_x < prepls_lims[1]))[0]
        pre, pre_x, header, pre_std = mufuns.get_asy(run, start_t = prepls_lims[0], end_t = prepls_lims[1], single_bin=True, year = run_year, alpha = run_alpha)


        strttime = 700
        intglen = 3000
        seq_offs = 0
        
        # get the post-pulse part
        postpls_lims = np.array([strttime, strttime+intglen]) + seq_offs
        postpls = np.where((asy_x > postpls_lims[0]) & (asy_x < postpls_lims[1]))[0]
        post, post_x, header, post_std = mufuns.get_asy(run, start_t = postpls_lims[0], end_t = postpls_lims[1], single_bin=True, year = run_year, alpha = run_alpha)

        
        
        amps[ii_run,0] = pre[0,asy_id]
        amps[ii_run,1] = first[0,asy_id]
        amps[ii_run,2] = post[0,asy_id]

        stds[ii_run,0] = pre_std[0,asy_id]
        stds[ii_run,1] = first_std[0,asy_id]
        stds[ii_run,2] = post_std[0,asy_id]
        
        # the center of the window. This is most relevant for interpls, where the time center is varying due to different window lengths!
        amps_x[ii_run,0] = pre_x
        amps_x[ii_run,1] = first_x
        amps_x[ii_run,2] = post_x
    
        
     

#%% plotting of integrals

# x axis: just use run number here
x_axis = exp_c[0]


plt.figure(26)
plt.clf()
plt.errorbar(x_axis, amps[:,1], yerr = stds[:,1], label='post')
plt.errorbar(x_axis, amps[:,0], yerr = stds[:,0], label= 'pre')
plt.errorbar(x_axis, amps[:,2], yerr = stds[:,2], label = 'post longer')

plt.legend(loc=0)
plt.xlabel('$t$ [ns]')
plt.ylabel('$A$')




plt.figure(27)
plt.clf()
plt.errorbar(x_axis, amps[:,1], yerr = stds[:,1], marker='o', linestyle='')

plt.xlabel('$t$ [ns]')
plt.ylabel('LF asymmetry')

# gsi = plt.gcf().get_size_inches()
# plt.gcf().set_size_inches(gsi[0]*0.75, gsi[1])

plt.gcf().set_size_inches(4.8, 3.2)
plt.rcParams["font.family"] = "sans-serif"


np.savez('S_260K_rotecho_1240G_toggle', pulspars=pulspars, amps = amps, stds = stds)
