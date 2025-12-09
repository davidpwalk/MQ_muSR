# -*- coding: utf-8 -*-
"""
Created on Tue May 27 15:10:31 2025

@author: doll_a
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants


# view those spectra



lt = np.load('../../Data/STO_17K_spectrum.npz')
ht = np.load('../../Data/STO_70K_spectrum.npz')


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

lt = np.load('../../Data/TS_15K_spectrum.npz')
ht = np.load('../../Data/GPS/TS_260K_spectrum.npz')

plt.show()
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



lt = np.load('../../Data/STO_17K_spectrum.npz')
ht = np.load('../../Data/GPS/TS_260K_spectrum.npz')

plt.show()
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
td = np.load('../../Data/STO_tdomain_dec2.npz')

plt.show()
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

#%% S-1 experiments, 2025, all at 260 K

TS1 = np.load('../../Data/GPS/TS_260K_spectrum.npz')
S1_setup = np.load('../../Data/S_260K_1stspec_nodegrader.npz')
S1 = np.load('../../Data/S_260K_spec.npz')

plt.show()
plt.figure(111)
plt.clf()
plt.errorbar(S1['B_vec'], S1['amps'][:,1], yerr = S1['stds'][:,1], marker='o', linestyle='', color='C0', label='S1' ,  solid_capstyle='projecting', capsize=3, markersize = 2.5)
plt.errorbar(TS1['B_vec'], TS1['amps'][:,1], yerr = TS1['stds'][:,1], marker='o', linestyle='', color='C1', label='TS1',  solid_capstyle='projecting', capsize=3, markersize = 2.5 )
plt.errorbar(S1_setup['B_vec'], S1_setup['amps'][:,1], yerr = S1_setup['stds'][:,1], marker='o', linestyle='', color='C2', label='S1 setup' ,  solid_capstyle='projecting', capsize=3, markersize = 2.5)


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

#%% pulse duration exp


tswp = np.load('../../Data/S_260K_rotecho_tswp.npz')


normpls_ids = np.where(tswp['pulspars'][:,3]==0)[0]
rotecho_ids = np.where(tswp['pulspars'][:,3]==1)[0]


normpls_ids = np.where((tswp['pulspars'][:,3]==0)&(tswp['pulspars'][:,4]==-25))[0]
rotecho_ids = np.where((tswp['pulspars'][:,3]==1)&(tswp['pulspars'][:,4]==-25))[0]


normpls13_ids = np.where((tswp['pulspars'][:,3]==0)&(tswp['pulspars'][:,4]==-13))[0]
rotecho13_ids = np.where((tswp['pulspars'][:,3]==1)&(tswp['pulspars'][:,4]==-13))[0]


# plt.figure(116)
# plt.clf()
# plt.errorbar(tswp['pulspars'][normpls_ids,2], tswp['amps'][normpls_ids,1], yerr = tswp['stds'][normpls_ids,1], label='post')
# plt.errorbar(tswp['pulspars'][normpls_ids,2], tswp['amps'][normpls_ids,0], yerr = tswp['stds'][normpls_ids,0], label= 'pre')
# plt.errorbar(tswp['pulspars'][normpls_ids,2], tswp['amps'][normpls_ids,2], yerr = tswp['stds'][normpls_ids,2], label = 'post longer')

# plt.legend(loc=0)
# plt.xlabel('$t$ [ns]')
# plt.ylabel('$A$')



plt.show()
plt.figure(117)
plt.clf()
plt.errorbar(tswp['pulspars'][normpls_ids,2], tswp['amps'][normpls_ids,1], yerr = tswp['stds'][normpls_ids,1], marker='o', linestyle='', color='C2', label = 'norm' )
plt.errorbar(tswp['pulspars'][rotecho_ids,2], tswp['amps'][rotecho_ids,1], yerr = tswp['stds'][rotecho_ids,1], marker='o', linestyle='', color='C1', label = 'rotecho' )


plt.legend(loc=0)
plt.xlabel('$t_\mathrm{pulse}$ [ns]')
plt.ylabel('LF asymmetry')


plt.show()
plt.figure(118)
plt.clf()
plt.errorbar(tswp['pulspars'][normpls13_ids,2], tswp['amps'][normpls13_ids,1], yerr = tswp['stds'][normpls13_ids,1], marker='o', linestyle='', color='C2', label = 'norm13' )
plt.errorbar(tswp['pulspars'][rotecho13_ids,2], tswp['amps'][rotecho13_ids,1], yerr = tswp['stds'][rotecho13_ids,1], marker='o', linestyle='', color='C1', label = 'rotecho13' )


plt.legend(loc=0)
plt.xlabel('$t_\mathrm{pulse}$ [ns]')
plt.ylabel('LF asymmetry')


plt.show()
plt.figure(119)
plt.clf()
plt.errorbar(tswp['pulspars'][normpls_ids,2], tswp['amps'][normpls_ids,1], yerr = tswp['stds'][normpls_ids,1], marker='o', linestyle='', color='C2', label = 'norm' )
plt.errorbar(tswp['pulspars'][rotecho_ids,2], tswp['amps'][rotecho_ids,1], yerr = tswp['stds'][rotecho_ids,1], marker='o', linestyle='', color='C1', label = 'rotecho' )

plt.errorbar(tswp['pulspars'][normpls13_ids,2], tswp['amps'][normpls13_ids,1], yerr = tswp['stds'][normpls13_ids,1], marker='o', linestyle='', color='C4', label = 'norm13' )
plt.errorbar(tswp['pulspars'][rotecho13_ids,2], tswp['amps'][rotecho13_ids,1], yerr = tswp['stds'][rotecho13_ids,1], marker='o', linestyle='', color='C3', label = 'rotecho13' )


plt.legend(loc=0)
plt.xlabel('$t_\mathrm{pulse}$ [ns]')
plt.ylabel('LF asymmetry')

#%% pulse amplitude sweeps


ampswp = np.load('../../Data/S_260K_rotecho_ampswp.npz')

normpls_ids = np.where(ampswp['pulspars'][:,3]==0)[0]
rotecho_ids = np.where(ampswp['pulspars'][:,3]==1)[0]


# plt.figure(156)
# plt.clf()
# plt.errorbar(ampswp['pulspars'][normpls_ids,4], ampswp['amps'][normpls_ids,1], yerr = ampswp['stds'][normpls_ids,1], label='post')
# plt.errorbar(ampswp['pulspars'][normpls_ids,4], ampswp['amps'][normpls_ids,0], yerr = ampswp['stds'][normpls_ids,0], label= 'pre')
# plt.errorbar(ampswp['pulspars'][normpls_ids,4], ampswp['amps'][normpls_ids,2], yerr = ampswp['stds'][normpls_ids,2], label = 'post longer')

# plt.legend(loc=0)
# plt.xlabel('Drive [dBm]')
# plt.ylabel('$A$')



plt.show()
plt.figure(157)
plt.clf()
plt.errorbar(ampswp['pulspars'][normpls_ids,4], ampswp['amps'][normpls_ids,1], yerr = ampswp['stds'][normpls_ids,1], marker='o', linestyle='', color='C2', label = 'norm' )
plt.errorbar(ampswp['pulspars'][rotecho_ids,4], ampswp['amps'][rotecho_ids,1], yerr = ampswp['stds'][rotecho_ids,1], marker='o', linestyle='', color='C1', label = 'rotecho' )

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


#%% On/Off to check satelite peak at 1240 G, versus 1250 G
onsat = np.load('../../Data/S_260K_rotecho_1240G_toggle.npz')
offsat = np.load('../../Data/S_260K_rotecho_1250G_toggle.npz')

# stitch together the data for now
amps = np.vstack((onsat['amps'],offsat['amps']))
stds = np.vstack((onsat['stds'],offsat['stds']))
x_axis =  np.vstack((onsat['runs'][:,None],offsat['runs'][:,None]))



# plt.figure(136)
# plt.clf()
# plt.errorbar(x_axis, amps[:,1], yerr = stds[:,1], label='post')
# plt.errorbar(x_axis, amps[:,0], yerr = stds[:,0], label= 'pre')
# plt.errorbar(x_axis, amps[:,2], yerr = stds[:,2], label = 'post longer')

# plt.legend(loc=0)
# plt.xlabel('run number')
# plt.ylabel('$A$')



plt.show()
plt.figure(137)
plt.clf()
plt.errorbar(x_axis, amps[:,1], yerr = stds[:,1], marker='o', linestyle='')
plt.errorbar(x_axis, amps[:,2], yerr = stds[:,2], marker='o', linestyle='')

plt.xlabel('run number')
plt.ylabel('LF asymmetry')

# gsi = plt.gcf().get_size_inches()
# plt.gcf().set_size_inches(gsi[0]*0.75, gsi[1])

plt.gcf().set_size_inches(4.8, 3.2)
plt.rcParams["font.family"] = "sans-serif"
plt.show()

