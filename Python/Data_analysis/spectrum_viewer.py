# -*- coding: utf-8 -*-
"""
Created on Tue May 27 15:10:31 2025

@author: doll_a
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants


# view those spectra
lt = np.load('STO_17K_spectrum.npz')
ht = np.load('STO_70K_spectrum.npz')


plt.figure(10)
plt.clf()
plt.errorbar(lt['B_vec'], lt['amps'][:,1], yerr = lt['stds'][:,1], marker='o', linestyle='', color='C0', label='17 K' )
plt.errorbar(ht['B_vec'], ht['amps'][:,1], yerr = ht['stds'][:,1], marker='o', linestyle='', color='C1', label='70 K' )

y_pos = 0.22
perc = 2

diff = y_pos/100*perc

y_bar = np.array([y_pos, y_pos-diff])
x_pos = np.array([1.0,1.0])*1195.0

plt.plot(x_pos, y_bar)

plt.legend(loc=0)
plt.xlabel('Magnetic Field [G]')
plt.ylabel('LF asymmetry')

# gsi = plt.gcf().get_size_inches()
# plt.gcf().set_size_inches(gsi[0]*0.75, gsi[1])

plt.gcf().set_size_inches(4.8, 3.2)
# plt.rcParams["font.family"] = "sans-serif"

plt.yticks([0.215, 0.22, 0.225, 0.23])


lt = np.load('TS_15K_spectrum.npz')
ht = np.load('TS_260K_spectrum.npz')


plt.figure(11)
plt.clf()
plt.errorbar(lt['B_vec'], lt['amps'][:,1], yerr = lt['stds'][:,1], marker='o', linestyle='', color='C0', label='15 K' ,  solid_capstyle='projecting', capsize=3, markersize = 2.5)
plt.errorbar(ht['B_vec'], ht['amps'][:,1], yerr = ht['stds'][:,1], marker='o', linestyle='', color='C1', label='260 K',  solid_capstyle='projecting', capsize=3, markersize = 2.5 )

y_pos = 0.21
perc = 2

diff = y_pos/100*perc

y_bar = np.array([y_pos, y_pos-diff])
x_pos = np.array([1.0,1.0])*1280

plt.plot(x_pos, y_bar)

# also add the linewidth

y_pos = 0.215
x_pos = 1235
x_w = 20

y_bar = np.array([y_pos, y_pos])
x_bar = np.array([x_pos-x_w/2, x_pos+x_w/2])

plt.plot(x_bar, y_bar)


plt.legend(loc=0)
plt.xlabel('Magnetic Field [G]')
plt.ylabel('LF asymmetry')

# gsi = plt.gcf().get_size_inches()
# plt.gcf().set_size_inches(gsi[0]*0.75, gsi[1])

plt.gcf().set_size_inches(4.8, 3.2)
# plt.rcParams["font.family"] = "sans-serif"

# sca_proposal = 0.6
# plt.gcf().set_size_inches(4.8 * sca_proposal, 3.2 * sca_proposal)
# plt.rcParams["font.family"] = "sans-serif"



lt = np.load('STO_17K_spectrum.npz')
ht = np.load('TS_260K_spectrum.npz')

plt.figure(12)
plt.clf()
plt.errorbar(lt['B_vec']-1182, lt['amps'][:,1], yerr = lt['stds'][:,1], marker='o', linestyle='', color='C2', label='15 K' ,  solid_capstyle='projecting', capsize=3, markersize = 2.5)
plt.errorbar(ht['B_vec']-1235, ht['amps'][:,1], yerr = ht['stds'][:,1], marker='o', linestyle='', color='C1', label='260 K',  solid_capstyle='projecting', capsize=3, markersize = 2.5 )

y_pos = 0.21
perc = 2

diff = y_pos/100*perc

y_bar = np.array([y_pos, y_pos-diff])
x_pos = np.array([1.0,1.0])*(1280-1235)

plt.plot(x_pos, y_bar)

# also add the linewidth

y_pos = 0.215
x_pos = 0
x_w = 20

y_bar = np.array([y_pos, y_pos])
x_bar = np.array([x_pos-x_w/2, x_pos+x_w/2])

# plt.plot(x_bar, y_bar)


# plt.legend(loc=0)
plt.xlabel('Field Offset [G]')
plt.ylabel('LF asymmetry')

plt.xlim(-50,50)




#%% get the g-values

# get gamma
hbar = constants.hbar
e = constants.elementary_charge
me = constants.electron_mass

mub = e*hbar/(2*me)

g = 2.0035
# g = 2.0038
# g = 2.00261/0.015


g = 2.00
# g = 1.93


gamma_radians = g * mub/hbar
gamma_Hz = gamma_radians/(2*np.pi)
gamma_MHz = gamma_Hz * 1e-6
gamma_MHzG = gamma_MHz * 1e-4


f_uw = 3456

print('g = {}:\t{} G for {} MHz'.format(g, f_uw/gamma_MHzG, f_uw))


#%% also add the timedomain on-off resonance data. otherwise, always re-running the full eval script takes forever...


xlims = [0, 2200]


# td = np.load('STO_tdomain.npz')
td = np.load('STO_tdomain_dec2.npz')

plt.figure(50); plt.clf()
plt.errorbar(td['tdx'], td['tdy_off'], yerr=td['tdy_off_std'], marker='o', color = 'gray', linestyle='None')
plt.errorbar(td['tdx'], td['tdy_on'], yerr=td['tdy_on_std'], marker='o', color = 'black', linestyle='None')





plt.xlabel('$t$ [ns]')
plt.ylabel('$A(t)$')
plt.xlim(xlims)
plt.legend(['off', 'on'])


ssi = np.array([4.8, 3.2])
plt.gcf().set_size_inches(ssi*0.6)
# plt.rcParams["font.family"] = "sans-serif"
