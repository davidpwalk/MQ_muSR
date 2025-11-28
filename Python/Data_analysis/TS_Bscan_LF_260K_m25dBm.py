# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 16:29:27 2020

@author: doll_a
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy
import Python.mufuns as mufuns
from scipy import interpolate

from scipy import signal
from scipy.fftpack import fft, fftshift, fftfreq
from scipy import optimize
import iminuit
from iminuit import cost


exp_c = []

# exp_c.append(np.arange(3867,3869)) # alpha runs

exp_c.append(np.arange(3955, 3979))


run_year = 2024
# alpha taken from run 2432
run_alpha = [1.175, 1.10, 1.10]


for ii_exp, runs in enumerate(exp_c):

    B_vec = np.zeros_like(runs, dtype=float)

    amps = np.zeros((len(runs),3))
    stds = np.zeros((len(runs),3))
    amps_x = np.zeros((len(runs),3))
    


    fdxs = []
    fdys = []

    tdxs = []
    tdys = []
    
    asys = []
    asys_x = []
    evranges = []
    
    asys_dec = []
    asys_dec_x = []
    asys_dec_std = []

    # not used, but append_new only works with lists...
    pulspars = np.zeros((len(runs),4))
    

    for ii_run, run in enumerate(runs):

        # get the data, but actually more for the purpose of the header
        asy, asy_x, header, asy_err = mufuns.get_asy(run,run_year)

        B_vec[ii_run] = float(header['field'][:-4])

        # try to fetch info about pulse delays
        title = header['title']

        asy_id = 0
        dec = 144
        
        # get the integral of the entire curve
        # integration window
        strttime = 1300
        intglen = 1000
        seq_offs = 50
        
        
        firstpls_lims = np.array([strttime, strttime+intglen]) + seq_offs
        firstpls = np.where((asy_x > firstpls_lims[0]) & (asy_x < firstpls_lims[1]))[0]
        first, first_x, header, first_std = mufuns.get_asy(run, start_t = firstpls_lims[0], end_t = firstpls_lims[1], single_bin=True, year = run_year, alpha = run_alpha)
        
        
        # get the pre-pulse part
        prepls_lims = [30, strttime-seq_offs]
        prepls = np.where((asy_x > prepls_lims[0]) & (asy_x < prepls_lims[1]))[0]
        pre, pre_x, header, pre_std = mufuns.get_asy(run, start_t = prepls_lims[0], end_t = prepls_lims[1], single_bin=True, year = run_year, alpha = run_alpha)


        
        # get the post-pulse part
        postpls_lims = [strttime+intglen + 250, 4000]
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
    
        ev_lims = [0, 4e3]
        # align such that the the decimated data have their datapoints coinciding with the pulse startpoint
        deltat = asy_x[1]-asy_x[0]
        dec_t = deltat*dec
        dec_t_half = dec_t/2 # to adjust for the center
        t_shift = np.mod(strttime + seq_offs-dec_t_half,dec_t)
        ev_lims[0] += t_shift
        
        
        asy_dec, asy_dec_x, header, asy_dec_std = mufuns.get_asy(run, start_t = ev_lims[0], end_t = ev_lims[1], binning=dec, year = run_year, alpha = run_alpha)



        # pls_lims = [750, 1750]
        pls_lims = firstpls_lims
        pls = np.where((asy_dec_x > pls_lims[0]) & (asy_dec_x < pls_lims[1]))[0]


        # select ev_ran and asymmetry
        tdx = asy_dec_x[pls]*1e-3 # usecs, MHz
        tdy = asy_dec[pls]

        # assure 2D array
        if len(tdy.shape) == 1:
            tdy_new = np.zeros((tdy.shape[0],1),dtype='complex_')
            tdy_new[:,0] = tdy
            tdy = tdy_new

        tdxs.append(tdx)
        tdys.append(tdy)
        
        asys.append(asy)
        asys_x.append(asy_x)
        evranges.append([prepls, firstpls])
        
        asys_dec.append(asy_dec)
        asys_dec_x.append(asy_dec_x)
        asys_dec_std.append(asy_dec_std)
        
     
#%% plotting of binned data

seltrcs = np.arange(len(asys_dec))


sca = np.mean(amps[:,0])
voffs = 0.1
# sca = 1

uniqueflds = np.unique(B_vec)
colors = plt.cm.RdBu_r(np.linspace(0,0.4,len(uniqueflds)))
colors = plt.cm.cool(np.linspace(0,1.0,len(uniqueflds)))
#colors = plt.cm.autumn(np.linspace(0,1.0,len(uniqueflds)))
#colors = plt.cm.hsv(np.linspace(0,1.0,len(uniqueflds)))

for ii_tr, seltrc in enumerate(seltrcs):  
    if ii_tr == 0:
        plt.figure(4)
        plt.clf()
        
    # plt.plot(asys_dec_x[seltrc], asys_dec[seltrc][:,asy_id])
    popt = '-'
    # plt.plot(asys_dec_x[seltrc],asys_dec[seltrc][:,asy_id]/sca+ B_vec[seltrc],popt, color = colors[np.where(B_vec[seltrc] == uniqueflds)[0][0]])
    plt.plot(asys_dec_x[seltrc],asys_dec[seltrc][:,asy_id]/sca + ii_tr * voffs,popt, color = colors[np.where(B_vec[seltrc] == uniqueflds)[0][0]])


plt.xlabel('$t$ [ns]')
plt.ylabel('$B_0$ [G]')
plt.xlim(100,2550)
B0span = 13
B0_cent = 823.75
# plt.ylim([B0_cent-B0span/2, B0_cent+B0span/2])

# simple explanation figure
seltrcs = [len(B_vec)-1,8]

cols = ['black','gray']  
cols = ['gray','black']  

for ii_tr, seltrc in enumerate(seltrcs):  
    if ii_tr == 0:
        plt.figure(400)
        plt.clf()
        
    # plt.plot(asys_dec_x[seltrc], asys_dec[seltrc][:,asy_id])
    popt = '-'
    plt.errorbar(asys_dec_x[seltrc],asys_dec[seltrc][:,asy_id],yerr=asys_dec_std[seltrc][:,asy_id], marker='o', color = cols[ii_tr])
    # plt.plot(asys_dec_x[seltrc],asys_dec[seltrc][:,asy_id]/sca + ii_tr * voffs,popt, color = colors[np.where(B_vec[seltrc] == uniqueflds)[0][0]])


plt.plot(np.array([286, 286]), np.array([0.1,0.14]))
plt.plot(np.array([286, 286])+seq_offs, np.array([0.1,0.14]))
plt.plot(np.array([286, 286])+seq_offs+1250, np.array([0.1,0.14]))
plt.xlabel('$t$ [ns]')
plt.ylabel('$A(t)$')
# plt.ylim(-0.001,0.145)
plt.xlim(0,2400)



#%% plotting of integrals
sort_ids = np.argsort(B_vec)


plt.figure(6)
plt.clf()
plt.errorbar(B_vec[sort_ids], amps[sort_ids,1], yerr = stds[sort_ids,1], label='pulse')
plt.errorbar(B_vec[sort_ids], amps[sort_ids,0], yerr = stds[sort_ids,0], label= 'pre')
plt.errorbar(B_vec[sort_ids], amps[sort_ids,2], yerr = stds[sort_ids,2], label = 'post')

plt.legend(loc=0)
plt.xlabel('$B_0$ [G]')
plt.ylabel('$A$')



# lt = np.load('TS_15K_spectrum.npz')
#
#
# plt.figure(7)
# plt.clf()
# plt.errorbar(B_vec[sort_ids], amps[sort_ids,1], yerr = stds[sort_ids,1], marker='o', linestyle='', color='C2' )
# plt.errorbar(lt['B_vec'], lt['amps'][:,1], yerr = lt['stds'][:,1], marker='o', linestyle='', color='C1' )
#
# y_pos = 0.22
# perc = 2
#
# diff = y_pos/100*perc
#
# y_bar = np.array([y_pos, y_pos-diff])
# x_pos = np.array([1.0,1.0])*1280.0
#
# plt.plot(x_pos, y_bar)
#
# plt.legend(loc=0)
# plt.xlabel('Magnetic Field [G]')
# plt.ylabel('LF asymmetry')
#
# # gsi = plt.gcf().get_size_inches()
# # plt.gcf().set_size_inches(gsi[0]*0.75, gsi[1])
#
# plt.gcf().set_size_inches(4.8, 3.2)
# plt.rcParams["font.family"] = "sans-serif"


np.savez('../../Data/TS_260K_spectrum', B_vec=B_vec[sort_ids], amps = amps[sort_ids,:], stds = stds[sort_ids,:])
