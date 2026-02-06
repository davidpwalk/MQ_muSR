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
3826 + 3827
exp_c.append(np.hstack((np.arange(3826, 3828),np.arange(3844,3864))))


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
            pulspars[ii_run,4] = -25.0 # default value from those where it wasn't put
            
            
            pstr = title.split(',')[5]
            pulspars[ii_run,4] = float(pstr.replace('m', '-').replace('dB',''))
            
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
normpls_ids = np.where(pulspars[:,3]==0)[0]
rotecho_ids = np.where(pulspars[:,3]==1)[0]


plt.figure(26)
plt.clf()
plt.errorbar(pulspars[normpls_ids,4], amps[normpls_ids,1], yerr = stds[normpls_ids,1], label='post')
plt.errorbar(pulspars[normpls_ids,4], amps[normpls_ids,0], yerr = stds[normpls_ids,0], label= 'pre')
plt.errorbar(pulspars[normpls_ids,4], amps[normpls_ids,2], yerr = stds[normpls_ids,2], label = 'post longer')

plt.legend(loc=0)
plt.xlabel('Drive [dBm]')
plt.ylabel('$A$')




plt.figure(27)
plt.clf()
plt.errorbar(pulspars[normpls_ids,4], amps[normpls_ids,1], yerr = stds[normpls_ids,1], marker='o', linestyle='', color='C2', label = 'norm' )
plt.errorbar(pulspars[rotecho_ids,4], amps[rotecho_ids,1], yerr = stds[rotecho_ids,1], marker='o', linestyle='', color='C1', label = 'rotecho' )

y_pos = 0.22
perc = 2

diff = y_pos/100*perc

y_bar = np.array([y_pos, y_pos-diff])
x_pos = np.array([1.0,1.0])*1280.0

# plt.plot(x_pos, y_bar)

plt.legend(loc=0)
plt.xlabel('Drive [dBm]')
plt.ylabel('LF asymmetry')

# gsi = plt.gcf().get_size_inches()
# plt.gcf().set_size_inches(gsi[0]*0.75, gsi[1])

plt.gcf().set_size_inches(4.8, 3.2)
plt.rcParams["font.family"] = "sans-serif"


np.savez('S_260K_rotecho_ampswp', pulspars=pulspars, amps = amps, stds = stds)