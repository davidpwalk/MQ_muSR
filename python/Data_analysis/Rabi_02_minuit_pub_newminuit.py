# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 16:29:27 2020

@author: doll_a
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy
import python.utils.mufuns as mufuns

from scipy.fftpack import fft, fftshift, fftfreq
import iminuit
from iminuit import cost


exp_c = []
# exp_c.append(np.arange(2154,2191))
exp_c.append(np.arange(2289,2303))


for ii_exp, runs in enumerate(exp_c):

    B_vec = np.zeros_like(runs)

    amps = np.zeros((len(runs),2))
    stds = np.zeros((len(runs),2))
    amps_x = np.zeros((len(runs),2))

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

    pulspars = np.zeros((len(runs),4))
    

    for ii_run, run in enumerate(runs):

        # get the data, but actually more for the purpose of the header
        asy, asy_x, header, asy_err = mufuns.get_asy(run)

        B_vec[ii_run] = float(header['field'][:-4])

        # try to fetch info about pulse delays
        title = header['title']

        asy_id = 0
        dec = 144
        
        # get the integral of the entire curve
        # integration window
        strttime = 286+5
        intglen = 1250
        seq_offs = 200
        
        
        firstpls_lims = np.array([strttime, strttime+intglen]) + seq_offs
        firstpls = np.where((asy_x > firstpls_lims[0]) & (asy_x < firstpls_lims[1]))[0]
        first, first_x, header, first_std = mufuns.get_asy(run, start_t = firstpls_lims[0], end_t = firstpls_lims[1], single_bin=True)
        
        
        # get the pre-pulse part
        prepls_lims = [100, strttime+seq_offs]
        prepls = np.where((asy_x > prepls_lims[0]) & (asy_x < prepls_lims[1]))[0]
        pre, pre_x, header, pre_std = mufuns.get_asy(run, start_t = prepls_lims[0], end_t = prepls_lims[1], single_bin=True)

        
        amps[ii_run,0] = pre[0,asy_id]
        amps[ii_run,1] = first[0,asy_id]

        stds[ii_run,0] = pre_std[0,asy_id]
        stds[ii_run,1] = first_std[0,asy_id]
        
        # the center of the window. This is most relevant for interpls, where the time center is varying due to different window lengths!
        amps_x[ii_run,0] = pre_x
        amps_x[ii_run,1] = first_x
    
        ev_lims = [0, 4e3]
        # align such that the the decimated data have their datapoints coinciding with the pulse startpoint
        deltat = asy_x[1]-asy_x[0]
        dec_t = deltat*dec
        dec_t_half = dec_t/2 # to adjust for the center
        t_shift = np.mod(strttime + seq_offs-dec_t_half,dec_t)
        ev_lims[0] += t_shift
        
        
        asy_dec, asy_dec_x, header, asy_dec_std = mufuns.get_asy(run, start_t = ev_lims[0], end_t = ev_lims[1], binning=dec)



        # pls_lims = [750, 1750]
        pls_lims = firstpls_lims
        pls = np.where((asy_dec_x > pls_lims[0]) & (asy_dec_x < pls_lims[1]))[0]


        # select ev_ran and asymmetry
        tdx = asy_dec_x[pls]*1e-3 # usecs, MHz
        tdy = asy_dec[pls]

        # assure 2D array
        if len(tdy.shape) == 1:
            tdy_new = np.zeros((tdy.shape[0], 1), dtype=np.complex128)
            tdy_new[:,0] = tdy
            tdy = tdy_new

        # fft
        DC_remove = True
        if DC_remove:
            tdy_w = tdy - np.mean(tdy,axis=0) # DC removal
        else:
            tdy_w = tdy


        N_fft = int(2**np.ceil(np.log2(2*len(tdy_w))))
        fdy = np.zeros((N_fft,tdy_w.shape[1]),dtype=np.complex128)
        for ii in np.arange(fdy.shape[1]):
            fdy[:,ii] = fftshift(fft(tdy_w[:,ii],axis=0,n=N_fft),axes=0) / tdy_w.shape[0]
        fdx = fftfreq(len(fdy),d=tdx[1]-tdx[0])
        fdx = fftshift(fdx)

        # correct phase, first order
        dfdx = fdx[1] - fdx[0]
        maid = abs(fdy).argmax(0)
        maid = int(np.mean(maid)) # take mean in case of multiple traces
        pharan_MHz = 2.0
        pharan_pts = int(np.round(pharan_MHz/2/dfdx))
        pharan = np.arange(maid - pharan_pts, maid + pharan_pts + 1)

        phase_zeroorder = True
        phase_firstorder = False
        if phase_zeroorder:
            for ii in np.arange(fdy.shape[1]):
                corrpha = np.angle(fdy[maid,ii])
                fdy[:,ii] = fdy[:,ii] * np.exp(-1j*corrpha)

                # TODO!
                if phase_firstorder:

                    phapol = np.polyfit(fdx[pharan],np.unwrap(np.angle(fdy[pharan,ii])),1)

                    totpha = np.polyval(phapol,fdx[maid])

                    tdy[:,ii] = tdy[:,ii]* np.exp(-1j*totpha)

                    tshift = -phapol[0] /2/np.pi

                    #print(tshift)
                    #fdy[:,ii] = fdy[:,ii] * np.exp(-1j*np.polyval(phapol,fdx))

        fdxs.append(fdx)
        fdys.append(fdy)

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
voffs = 0.4
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
    plt.plot(asys_dec_x[seltrc],asys_dec[seltrc][:,asy_id]/sca+ B_vec[seltrc],popt, color = colors[np.where(B_vec[seltrc] == uniqueflds)[0][0]])
    # plt.plot(asys_dec_x[seltrc],asys_dec[seltrc][:,asy_id]/sca + ii_tr * voffs,popt, color = colors[np.where(B_vec[seltrc] == uniqueflds)[0][0]])


plt.xlabel('$t$ [ns]')
plt.ylabel('$B_0$ [G]')
plt.xlim(100,2550)
B0span = 13
B0_cent = 823.75
plt.ylim([B0_cent-B0span/2, B0_cent+B0span/2])
plt.show()

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
plt.ylim(-0.001,0.145)
plt.xlim(0,2400)
plt.show()


# simple explanation figure
seltrcs = [len(B_vec)-1,8]

seltrcs = [8]
for ii_tr, seltrc in enumerate(seltrcs):  
    if ii_tr == 0:
        plt.figure(401)
        plt.clf()
        
    # plt.plot(asys_dec_x[seltrc], asys_dec[seltrc][:,asy_id])
    popt = '-'
    plt.errorbar(asys_dec_x[seltrc],asys_dec[seltrc][:,asy_id],yerr=asys_dec_std[seltrc][:,asy_id], marker='o')
    # plt.plot(asys_dec_x[seltrc],asys_dec[seltrc][:,asy_id]/sca + ii_tr * voffs,popt, color = colors[np.where(B_vec[seltrc] == uniqueflds)[0][0]])


# plt.plot(np.array([286, 286]), np.array([0.1,0.14]))
# plt.plot(np.array([286, 286])+seq_offs, np.array([0.1,0.14]))
# plt.plot(np.array([286, 286])+seq_offs+1250, np.array([0.1,0.14]))
plt.xlabel('$t$ [ns]')
plt.ylabel('$A(t)$')
plt.ylim(-0.001,0.145)
plt.show()


#%% FFT plot

sca = 0
for ii in range(len(fdxs)):
    sca = np.max([sca,abs(fdys[ii]).max()])
voffs = 0.2

uniqueflds = np.unique(B_vec)
colors = plt.cm.RdBu_r(np.linspace(0,0.4,len(uniqueflds)))
colors = plt.cm.cool(np.linspace(0,1.0,len(uniqueflds)))
#colors = plt.cm.autumn(np.linspace(0,1.0,len(uniqueflds)))
#colors = plt.cm.hsv(np.linspace(0,1.0,len(uniqueflds)))


pcol_dta = []

plt.figure(5)
plt.clf()
for ii in range(len(fdxs)):
    popt = '-'
    plt.plot(fdxs[ii],abs(fdys[ii][:,asy_id])/sca+ B_vec[ii],popt, color = colors[np.where(B_vec[ii] == uniqueflds)[0][0]])
    # plt.plot(fdxs[ii],abs(fdys[ii][:,1])/sca + ii*voffs,popt, color = colors[np.where(B_vec[ii] == uniqueflds)[0][0]])
    
#plt.legend(loc=0)

plt.xlabel('$\\nu_\mathrm{eff}$ [MHz]')
plt.ylabel('$B_0$ [G]')

plt.xlim([1,17])
plt.xlim([4,12])
B0span = 13
B0_cent = 823.75-0.5
plt.ylim([B0_cent-B0span/2,B0_cent+B0span/2])
plt.show()

#%% plotting of integrals
sort_ids = np.argsort(B_vec)


plt.figure(6)
plt.clf()
plt.errorbar(B_vec[sort_ids], amps[sort_ids,1], yerr = stds[sort_ids,1])
plt.errorbar(B_vec[sort_ids], amps[sort_ids,0], yerr = stds[sort_ids,0])

    
plt.xlabel('$B_0$ [G]')
plt.ylabel('$A$')
plt.show()

# store back the ratio
A_bg = amps[sort_ids,1] / amps[sort_ids,0]
A_bg_std = A_bg * np.sqrt(stds[sort_ids,1]**2/amps[sort_ids,1]**2 + stds[sort_ids,0]**2/amps[sort_ids,0]**2)


#%% get the traces normalized

trcs = []
trcs_std = []
trcs_x = []


# reobtain the pulse limits (in principle done already for tdx/tdy above, but here with the ratio)
firstpls_lims_dec = np.array([strttime, strttime+intglen]) + seq_offs

for ii in range(len(asys_dec)):  
    
    pls = np.where((asys_dec_x[ii] >= firstpls_lims_dec[0]-deltat/2) & (asys_dec_x[ii] < firstpls_lims_dec[1]))[0]
    trc = asys_dec[ii][pls,asy_id] / amps[ii,0]
    trc_x = asys_dec_x[ii][pls]
    trc_x = trc_x - trc_x[0]
    # stddev of ratio
    trc_std = trc*np.sqrt(asys_dec_std[ii][pls,asy_id]**2/asys_dec[ii][pls,asy_id]**2 + stds[ii,0]**2/amps[ii,0]**2)

    # trc = trc - np.mean(trc)

    trcs.append(trc)
    trcs_x.append(trc_x)
    trcs_std.append(trc_std)



uniqueflds = np.unique(B_vec)
colors = plt.cm.RdBu_r(np.linspace(0,0.4,len(uniqueflds)))
colors = plt.cm.cool(np.linspace(0,1.0,len(uniqueflds)))
#colors = plt.cm.autumn(np.linspace(0,1.0,len(uniqueflds)))
#colors = plt.cm.hsv(np.linspace(0,1.0,len(uniqueflds)))


offset = 0.2
for ii in range(len(trcs)):  
    if ii == 0:
        plt.figure(10)
        plt.clf()
        
    
    plt.errorbar(trcs_x[ii],trcs[ii] + ii*offset, yerr=trcs_std[ii], color = colors[np.where(B_vec[ii] == uniqueflds)[0][0]])


plt.xlabel('$t$ [ns]')
plt.ylabel('$A/A_0$')
plt.show()


#%% fit all of them

opts = []
fittrcs = []
fittrcs_fine = []

for ii in range(len(trcs)):  
    
    trc = trcs[ii]
    trc_x = trcs_x[ii]
    trc_err = trcs_std[ii]
    
    # fit function
    def fitfunction(x, trc_x = trc_x):
        fitfun = x[0] + x[1]*np.cos(2*np.pi*x[2]*trc_x + x[3])*np.exp(-trc_x/x[4])
        return fitfun
    
    def fitfunction_names(trc_x, A, Ap, nu, phi, tau):
        fitfun = A + Ap*np.cos(2*np.pi*nu*trc_x + phi)*np.exp(-trc_x/tau)
        return fitfun
    
    def fitfunction_minuit(x):
        fitfun = x[0] + x[1]*np.cos(2*np.pi*x[2]*trc_x + x[3])*np.exp(-trc_x/x[4])
        return fitfun
    
    fit_ran = np.where(trc_x < 11900)[0]
    
    # fit damped oscillation
    def damped_osc(x):
        return np.sum(abs(trc[fit_ran]-fitfunction(x)[fit_ran]))
    
    
    
    ang_ran = 20/180.0*np.pi
    ranges = [(-0.01, 0.01), (0.03, 0.45), (6e-3, 12e-3), (-ang_ran, +ang_ran), (0.2e3, 3e3)]
    ranges = [(-0.01, 1.5), (0.02, 0.6), (6e-3, 12e-3), (-ang_ran, +ang_ran), (0.2e3, 5e3)]
    bruteopt_N = 10
    bruteopt = scipy.optimize.brute(damped_osc, ranges, Ns = bruteopt_N)
    
    # minuit Xi squared initialization
    # costf = iminuit.cost.LeastSquares(trc_x, trc, trc_err, fitfunction)
    costf = iminuit.cost.LeastSquares(trc_x, trc, trc_err, fitfunction_names)
    m = iminuit.Minuit(costf, A = bruteopt[0], Ap = bruteopt[1], nu = bruteopt[2], phi = bruteopt[3], tau = bruteopt[4])  # minuit obj with initial value at bruteopt
    
    # costf = iminuit.cost.LeastSquares(trc_x, trc, trc_err, fitfunction_minuit)
    # m = iminuit.Minuit(costf, bruteopt)  # minuit obj with initial value at bruteopt
    
    
    m.migrad()
    m.hesse()
    m.minos()
    
    x0 = bruteopt
    bnds = ranges
    # opt = scipy.optimize.minimize(damped_osc, x0, method= 'L-BFGS-B', bounds = bnds)
    # opt = iminuit.minimize(damped_osc, x0, bounds = bnds)
    opt = iminuit.minimize(damped_osc, x0)
    
    vals = np.zeros((len(m.values),1))
    for ii in range(len(vals)):
        vals[ii] = m.values[ii]
        
    opt.x = vals
    opt.minuit = m
    
    
    
    # re-run until a valid minimum is found
    for ii in range(10):
        if opt.minuit.fmin.is_valid:
            print('valid minimum after {} additional runs'.format(ii))
            break
        opt.minuit.migrad()
    
    opt.minuit.hesse()
    opt.minuit.minos()
    
    # minimization
    #ranges = [(0.1e-12,0.1e-9), (0.1e-12,0.1e-9)]
    #x0 = [1e-12, 1e-12]
    x0 = opt.x
    
    minimizer_kwargs = {'method': 'L-BFGS-B', 'bounds': ranges}
    
    # function that checks whether random hop is within range 
    def bounds_for_randomhop(**kwargs):
         x = kwargs["x_new"]
         inrange = True
         for ii in np.arange(len(ranges)):
             tsts =  (ranges[ii][0] <= x[ii] <= ranges[ii][1])
             if not tsts:
                 print('jumped too far!')
             inrange = inrange and (ranges[ii][0] <= x[ii] <= ranges[ii][1])
         return inrange
    
    stepsize = 0.1e-12
    
    
    # glopt = scipy.optimize.basinhopping(damped_osc, x0, accept_test = bounds_for_randomhop, minimizer_kwargs=minimizer_kwargs, stepsize=stepsize, niter=200)
    glopt = opt
    
    opts.append([bruteopt, opt, glopt])
    
    trc_x_fine = np.linspace(trc_x[0], trc_x[-1], 1001)

    globfit_fine = fitfunction(glopt.x, trc_x = trc_x_fine)
    globfit = fitfunction(glopt.x)
    locfit = fitfunction(opt.x)
    locfit_fine = fitfunction(opt.x, trc_x = trc_x_fine)
    brutefit = fitfunction(bruteopt)
    brutefit_fine = fitfunction(bruteopt, trc_x = trc_x_fine)
    
    fittrcs.append([brutefit, locfit, globfit])
    fittrcs_fine.append([brutefit_fine, locfit_fine, globfit_fine, trc_x_fine])
    
#%% inspect results

x_brute = np.zeros((len(opts[0][0]),len(opts)))
x_opt = np.zeros_like(x_brute)
x_glo = np.zeros_like(x_brute)

x_opt_err = np.zeros((len(opts[0][0]),len(opts),2))
x_opt_errsrc = np.zeros_like(x_brute)

# extract the fitting parameters
for ii, sol in enumerate(opts):
    x_brute[:,ii] = sol[0]
    x_opt[:,ii] = sol[1].x[:,0]
    x_glo[:,ii] = sol[2].x[:,0]
    
# get the fitting errors
for ii, sol in enumerate(opts):
    imin = sol[1].minuit
    
    # populate with errors from Hessian
    hessret = imin.params
    for jj in range(len(hessret)):
        x_opt_err[jj, ii, 0] = hessret[jj].error
        x_opt_err[jj, ii, 1] = hessret[jj].error
    
    # print(x_opt_err[0,ii,:])
    # try to put minos errors
    try:
        minosret = imin.merrors
        for jj in range(len(minosret)):
            if minosret[jj]['lower_valid']:
                x_opt_err[jj, ii, 0] = -minosret[jj]['lower']
                x_opt_errsrc[jj,ii] += 1
               
            if minosret[jj]['upper_valid']:
                x_opt_err[jj, ii, 1] = minosret[jj]['upper']
                x_opt_errsrc[jj,ii] += 2
    except:
        pass
    
    # print(x_opt_err[jj,ii,:])
    

pids = np.arange(len(B_vec)-1)

plt.figure(19)
plt.clf()
# plt.plot(B_vec[pids], x_brute[0,pids])
# plt.plot(B_vec[pids], x_opt[0,pids])
plt.errorbar(B_vec[pids]/10, x_opt[0,pids], yerr = x_opt_err[0,pids,:].transpose(),color='black')
# plt.plot(B_vec[pids], x_glo[0,pids])
# plt.errorbar(B_vec[pids]/10, A_bg[pids], yerr=A_bg_std[pids], color='gray')
plt.xlabel('$B_0$ [mT]')
plt.ylabel('$A_\mathrm{bg}$')
plt.show()

plt.figure(20)
plt.clf()
# plt.plot(B_vec[pids], x_brute[1,pids])
plt.errorbar(B_vec[pids], x_opt[1,pids], yerr = x_opt_err[1,pids,:].transpose())
# plt.plot(B_vec[pids], x_opt[1,pids])
plt.plot(B_vec[pids], x_glo[1,pids])
plt.xlabel('$B_0$ [G]')
plt.ylabel('$\\Delta$$A$')
plt.show()

plt.figure(21)
plt.clf()
# plt.plot(B_vec[pids], x_brute[2,pids]*1e3)
# plt.plot(B_vec[pids], x_opt[2,pids]*1e3)
plt.errorbar(B_vec[pids], x_opt[2,pids]*1e3, yerr = x_opt_err[2,pids,:].transpose()*1e3)
plt.plot(B_vec[pids], x_glo[2,pids]*1e3)
plt.xlabel('$B_0$ [G]')
plt.ylabel('$\\nu_\mathrm{eff}$ [MHz]')
plt.show()

plt.figure(22)
plt.clf()
# plt.plot(B_vec[pids], x_brute[3,pids]/np.pi*180)
# plt.plot(B_vec[pids], x_opt[3,pids]/np.pi*180)
plt.errorbar(B_vec[pids], x_opt[3,pids]/np.pi*180, yerr = x_opt_err[3,pids,:].transpose()/np.pi*180)
plt.plot(B_vec[pids], x_glo[3,pids]/np.pi*180)
plt.xlabel('$B_0$ [G]')
plt.show()
    

plt.figure(23)
plt.clf()
# plt.plot(B_vec[pids], x_brute[4,pids])
# plt.plot(B_vec[pids], x_opt[4,pids])
plt.errorbar(B_vec[pids]/10, x_opt[4,pids], yerr = x_opt_err[4,pids,:].transpose(), color='black')
# plt.plot(B_vec[pids]/10, x_glo[4,pids])
plt.xlabel('$B_0$ [mT]')
plt.ylabel('$\\tau_\mathrm{dec}$')
plt.show()

#%% try to fit weff to the fitted frequency

trc = x_opt[2,pids]*1e3
trc_x = B_vec[pids]
trc_std = x_opt_err[2,pids]*1e3

# fit function
def fitfunction2(x, trc_x = trc_x):
    fitfun = np.sqrt(x[0]**2 + (x[1]*(trc_x - x[2]))**2)
    return fitfun


def fitfunction_names2(trc_x, nu1, gamma, Bres):
    fitfun = np.sqrt(nu1**2 + (gamma*(trc_x - Bres))**2)
    return fitfun    

fit_ran = np.arange(len(trc))

# fit weff
def fit_weff(x):
    return np.sum(abs(trc[fit_ran]-fitfunction2(x)[fit_ran]))

ranges = [(6.0, 8.0), (0.0, 3.0), (823, 827)]
bruteopt_N = 20
bruteopt2 = scipy.optimize.brute(fit_weff, ranges, Ns = bruteopt_N)

costf2 = iminuit.cost.LeastSquares(trc_x, trc, np.max(trc_std,axis=1), fitfunction_names2)
m2 = iminuit.Minuit(costf2, nu1 = bruteopt2[0], gamma = bruteopt2[1], Bres = bruteopt2[2])  # minuit obj with initial value at bruteopt

# costf = iminuit.cost.LeastSquares(trc_x, trc, trc_err, fitfunction_minuit)
# m = iminuit.Minuit(costf, bruteopt)  # minuit obj with initial value at bruteopt


m2.migrad()
m2.hesse()
m2.minos()

x0 = bruteopt
bnds = ranges
# opt = scipy.optimize.minimize(damped_osc, x0, method= 'L-BFGS-B', bounds = bnds)
# opt = iminuit.minimize(damped_osc, x0, bounds = bnds)
opt2 = iminuit.minimize(fit_weff, x0)

vals = np.zeros((len(m2.values),1))
for ii in range(len(vals)):
    vals[ii] = m2.values[ii]
    
opt2.x = vals
opt2.minuit = m2

# minimization
#ranges = [(0.1e-12,0.1e-9), (0.1e-12,0.1e-9)]
#x0 = [1e-12, 1e-12]
x0 = opt2.x

minimizer_kwargs = {'method': 'L-BFGS-B', 'bounds': ranges}

# function that checks whether random hop is within range 
def bounds_for_randomhop(**kwargs):
     x = kwargs["x_new"]
     inrange = True
     for ii in np.arange(len(ranges)):
         tsts =  (ranges[ii][0] <= x[ii] <= ranges[ii][1])
         if not tsts:
             print('jumped too far!')
         inrange = inrange and (ranges[ii][0] <= x[ii] <= ranges[ii][1])
     return inrange

stepsize = 0.1e-12


glopt2 = scipy.optimize.basinhopping(fit_weff, x0[:,0], accept_test = bounds_for_randomhop, minimizer_kwargs=minimizer_kwargs,
stepsize=stepsize, niter=200)



trc_x_fine = np.linspace(trc_x[0], trc_x[-1], 1001)

globfit_fine = fitfunction2(glopt2.x, trc_x = trc_x_fine)
globfit = fitfunction2(glopt2.x)
locfit = fitfunction2(opt2.x)
locfit_fine = fitfunction2(opt2.x, trc_x = trc_x_fine)
brutefit = fitfunction2(bruteopt2)
brutefit_fine = fitfunction2(bruteopt2, trc_x = trc_x_fine)

# get the fitting errors
x_opt2_err = np.zeros((len(opt2.x),2))
x_opt2_errsrc = np.zeros_like(bruteopt2)


imin = m2
 
# populate with errors from Hessian
hessret = imin.params
for jj in range(len(hessret)):
    x_opt2_err[jj, 0] = hessret[jj].error
    x_opt2_err[jj, 1] = hessret[jj].error
 
# print(x_opt_err[0,ii,:])
# try to put minos errors
try:
    minosret = imin.merrors
    for jj in range(len(minosret)):
        if minosret[jj].lower_valid:
            x_opt2_err[jj, 0] = -minosret[jj].lower
            x_opt2_errsrc[jj] += 1
        
        if minosret[jj].upper_valid:
            x_opt2_err[jj, 1] = minosret[jj].upper
            x_opt2_errsrc[jj] += 2
except:
    pass
 
 # print(x_opt_err[jj,ii,:])
 
# nu1, gamma, Bres
#     fitfun = np.sqrt(nu1**2 + (gamma*(trc_x - Bres))**2)
# get the detuning:
omega_fine = np.sqrt(opt2.x[1]*(trc_x_fine - opt2.x[2])**2)
omega = np.sqrt(opt2.x[1]*(trc_x - opt2.x[2])**2)

dA_mod = (1-(omega/locfit)**2)

dA_ampvec = np.linspace(0.35,0.45,101)
dA_ampvec = np.linspace(0.372,0.380,101)

dA_rmsd = np.zeros_like(dA_ampvec)

for ii in range(len(dA_rmsd)):
    dA_rmsd[ii] = np.sqrt(np.sum((dA_ampvec[ii]* dA_mod - x_opt[1,pids])**2))

rms_dA = dA_ampvec[np.argmin(dA_rmsd)]

# plt.figure(221)
# plt.clf()
# plt.plot(trc_x_fine/10.0, omega)

plt.figure(222)
plt.clf()
plt.plot(dA_ampvec, dA_rmsd)
plt.show()

plt.figure(25)
plt.clf()
plt.errorbar(B_vec[pids]/10.0, x_opt[2,pids]*1e3, yerr = x_opt_err[2,pids,:].transpose()*1e3,marker='o',linestyle='None', color='black')
plt.plot(trc_x_fine/10.0, locfit_fine,color='gray')
plt.xlabel('$B_0$ [mT]')
plt.ylabel('$\\nu_\mathrm{eff}$ [MHz]')
plt.show()

# plot the amplitude with the fitted scaling
plt.figure(26)
plt.clf()
# plt.plot(trc_x_fine/10,0.35*opt2.x[0]/locfit_fine,color='red')
plt.plot(trc_x_fine/10,rms_dA*(1-(omega_fine/locfit_fine)**2),color='gray')
# plt.plot(B_vec[pids], x_opt[1,pids])
plt.errorbar(B_vec[pids]/10, x_opt[1,pids], yerr = x_opt_err[1,pids,:].transpose(), color='black', linestyle = 'None', marker='o',)
plt.xlabel('$B_0$ [mT]')
plt.ylabel('$\\Delta$$A$')
plt.show()

# also plot Omega
def omega_fun(x, trc_x = trc_x):
    omega = (x[1]*(trc_x - x[2]))
    return omega

plt.figure(27)
plt.clf()
plt.plot(trc_x_fine, omega_fun(glopt2.x, trc_x = trc_x_fine))
omega_822 = omega_fun(glopt2.x, trc_x = 822)
plt.show()

# same for the constant background, which does depend on dA
psum_vec = np.linspace(0.7,1.2,101)
psum_vec = np.linspace(0.935,0.94,101)

psum_rmsd = np.zeros_like(psum_vec)

for ii in range(len(psum_rmsd)):
    psum_rmsd[ii] = np.sqrt(np.sum(((psum_vec[ii] - rms_dA*(1-(omega/locfit)**2))- x_opt[0,pids])**2))

rms_psum = psum_vec[np.argmin(psum_rmsd)]


plt.figure(223)
plt.clf()
plt.plot(psum_vec, psum_rmsd)
plt.show()

plt.figure(28)
plt.clf()
# plt.plot(B_vec[pids], x_brute[0,pids])
# plt.plot(B_vec[pids], x_opt[0,pids])
plt.plot(trc_x_fine/10, rms_psum - (rms_dA*(1-(omega_fine/locfit_fine)**2)),color='gray')
plt.errorbar(B_vec[pids]/10, x_opt[0,pids], yerr = x_opt_err[0,pids,:].transpose(),color='black', linestyle = 'None', marker='o')
# plt.plot(B_vec[pids], x_glo[0,pids])
# plt.errorbar(B_vec[pids]/10, A_bg[pids], yerr=A_bg_std[pids], color='gray')
plt.xlabel('$B_0$ [mT]')
plt.ylabel('$A_\mathrm{bg}$')
plt.show()

#%% Plot fitted Rabi curves

ii_sel = np.where(B_vec == 825)[0][0]
# ii_sel = np.where(B_vec == 826)[0][0]
# ii_sel = np.where(B_vec == 827)[0][0]
seltrcs = np.arange(5,12)

seltrcs = np.arange(len(B_vec))

seltrcs = np.flipud(seltrcs)

colors = plt.cm.cool(np.linspace(0,1.0,len(seltrcs)))


plt.figure(30)
plt.clf()   
offset = -0.5

for ii_off, ii_sel in enumerate(seltrcs):
    
    plt.errorbar(trcs_x[ii_sel],trcs[ii_sel] + ii_off*offset, yerr=trcs_std[ii_sel], color = colors[ii_off], linestyle = 'None', marker='o', label='{:.3g}'.format(B_vec[ii_sel]*1e-1))
    plt.plot(fittrcs_fine[ii_sel][3], fittrcs_fine[ii_sel][2] + ii_off*offset, color = colors[ii_off], linestyle=':', label = '_')
    plt.xlabel('$t$ [ns]')
    plt.ylabel('$A/A_0$')

plt.legend(loc=0)
plt.show()


# for the SI fig with all the traces
# plt.style.use('paper1x2h')
# then
plt.gcf().set_size_inches([6.4 , 8.61])

    
# np.savez('Rabi02_dta', trcs_x = trcs_x, trcs = trcs, trcs_std = trcs_std, B_vec = B_vec, fittrcs_fine=fittrcs_fine)



seltrcs = [5]

colors = ['black']


plt.figure(31)
plt.clf()   
offset = -0.5

for ii_off, ii_sel in enumerate(seltrcs):
    
    plt.errorbar(trcs_x[ii_sel],trcs[ii_sel] + ii_off*offset, yerr=trcs_std[ii_sel], color = colors[ii_off], linestyle = 'None', marker='o', label='{:.3g}'.format(B_vec[ii_sel]*1e-1))
    plt.plot(fittrcs_fine[ii_sel][3], fittrcs_fine[ii_sel][2] + ii_off*offset, color = colors[ii_off], linestyle=':', label = '_')
    plt.xlabel('$t$ [ns]')
    plt.ylabel('$A/A_0$')

plt.legend(loc=0)
plt.show()


    
#%% get omega at 822 G, to compare to Ramsey

ii_sel = np.where(B_vec == 827)[0][0]

B_ram = 822;

W_ram = 2.472
W_ram_e = 0.006

gamma = opt2.x[1]
Bres = opt2.x[2]

gamma_e = x_opt2_err[1,:].max()
Bres_e = x_opt2_err[2,:].max()

W_rabi = gamma * (-B_ram + Bres)

# error propagation assuming statistical independence
V_gamma = gamma_e**2
V_Bres = Bres_e**2

# re-parametrize: dB = B_ram - Bres
u_dB = B_ram - Bres
V_dB = Bres_e**2

Var_W =  u_dB**2 * V_gamma + gamma**2 * V_dB + V_dB*V_gamma
std_W = np.sqrt(Var_W) 

print('Omega at {} G from Rabi is {:.03f} +- {:.03f}'.format(B_ram, W_rabi[0], std_W[0]))
print('Omega at {} G from Ramsey is {} +- {}'.format(B_ram, W_ram, W_ram_e))


print('Bres is {:.03f} +- {:.03f} mT'.format(Bres[0]*1e-1, Bres_e*1e-1))

# also add the other relevant parameters, since it is annoying to dig around to find them

print('nu1 is {:.02f} +- {:.02f} MHz'.format(opt2.x[0][0], x_opt2_err[0,:].max()))
print('gamma^hat is {:.02f} +- {:.02f} MHz/T'.format(opt2.x[1][0]*1e1, x_opt2_err[1,:].max()*1e1)) # MHz/G to MHz/mT

# for B1, we need gamma_23 at the particular field value of 825 G from the SiO2_theory script: 7.235 MHz/mT
gamma_23 = 7.235
print('B_1 is {:.03f} +- {:.03f} MHz'.format(opt2.x[0][0]/gamma_23, x_opt2_err[0,:].max()/gamma_23))



#%% re-iterate on the very last trace, which has no rabis --> Fit T1 decay
ii = len(trcs)-1
trc = trcs[ii]
trc_x = trcs_x[ii]
trc_err = trcs_std[ii]

trc = asys_dec[ii][:,asy_id]
trc_x = asys_dec_x[ii]
trc_err = asys_dec_std[ii][:,asy_id]

    # plt.errorbar(asys_dec_x[seltrc],asys_dec[seltrc][:,asy_id],yerr=asys_dec_std[seltrc][:,asy_id], marker='o')



# fit function

def fitfunction_T1(trc_x, Ap, tau):
    fitfun = Ap*np.exp(-trc_x/tau)
    return fitfun

# minuit Xi squared initialization
costf = iminuit.cost.LeastSquares(trc_x, trc, trc_err, fitfunction_T1)
m = iminuit.Minuit(costf, Ap = 0.3, tau = 20e3)  # minuit obj with initial value at bruteopt

# costf = iminuit.cost.LeastSquares(trc_x, trc, trc_err, fitfunction_minuit)
# m = iminuit.Minuit(costf, bruteopt)  # minuit obj with initial value at bruteopt


m.migrad()
m.hesse()
m.minos()

#--> gives T1 on the order of 20 us

trc_fit = fitfunction_T1(trc_x,m.values['Ap'],m.values['tau'])

# plot
plt.figure(2000)
plt.clf()
plt.errorbar(trc_x, trc , yerr=trc_err)
plt.plot(trc_x, trc_fit)
plt.ylim(-0.1,1.1)
plt.show()


# also add to figure (400)

plt.figure(400)
plt.plot(trc_x, trc_fit)
plt.show()