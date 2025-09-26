# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 09:02:56 2020

@author: doll_a
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from scipy import linalg

from Python import jacobi_math

#%% analytical
B_vec = np.linspace(0.04,0.068,100)
B_vec = np.linspace(0.0,12.0,300)
# B_vec = np.linspace(0.0,0.9,300)
# B_vec = np.linspace(0.009,0.011,100)
# B_vec = np.linspace(0.04,0.09,200)

gn = 0.04257748
ge = 28.02495
w1 = ge*B_vec
w2 = -0.1355*B_vec
J = 2.008
#J = 4.4956 - 400e-3


w_mid = (w1 + w2)/2.0
w_diff = (w1 - w2)/2.0


alpha = 0.5*np.arctan2(J,w1-w2)
S =  np.sqrt((w_diff)**2.0 + (J/2)**2.0)
w12 = w_mid + J/2 - S
w13 = w_mid + J/2 + S
w24 = w_mid - J/2 + S
w34 = w_mid - J/2 - S

i12 = (1-np.sin(2*alpha))/4
i13 = (1+np.sin(2*alpha))/4
i24 = (1+np.sin(2*alpha))/4
i34 = (1-np.sin(2*alpha))/4


E1 = w_mid + J/2
E2 = -J/2 + S
E3 = -J/2 - S
E4 = -w_mid + J/2

w23 = E2 - E3


plt.figure(100)
plt.clf()
#plt.plot(B_vec, abs(w12))
#plt.plot(B_vec, abs(w13))
#plt.plot(B_vec, abs(w24))
plt.plot(B_vec, abs(w34))
#plt.plot(B_vec, abs(w23))

B_res_idx = np.where(abs(w34)<3.629)[0][0]
print(B_vec[B_res_idx])



plt.figure(101)
plt.clf()
plt.subplot(311)
plt.plot(B_vec, E1, label = '1')
plt.plot(B_vec, E2, label = '2')
plt.plot(B_vec, E3, label = '3')
plt.plot(B_vec, E4, label = '4')
plt.legend(loc=0)

plt.subplot(312)
plt.plot(B_vec, w12, label = '12')
plt.plot(B_vec, w13, label = '13')
plt.plot(B_vec, w24, label = '24')
plt.plot(B_vec, abs(w34), label = '34')
plt.legend(loc=0)

plt.subplot(313)
plt.plot(B_vec, i12, label = '12')
plt.plot(B_vec, i13, label = '13')
plt.plot(B_vec, i24, label = '24')
plt.plot(B_vec, i34, label = '34')
plt.legend(loc=0)

#%% give it a try with qutip

import qutip as qt
import Python.jacobi_math

B_vec = np.linspace(0.04,0.09,200)
# B_vec = np.array([0.0610])
B_vec = np.array([0.0822])
B_vec = np.array([0.188])
# B_vec = np.array([0.05])
# B_vec = np.array([0.003])
# B_vec = np.array([0])

# B_vec = np.linspace(0.0,0.1,300)
B_vec = np.linspace(0.0,0.9,300)
B_vec = np.linspace(0.0,0.6,300)
B_vec = np.linspace(0.0,0.6,51)
B_vec = np.linspace(0.0,0.6,151) # used for the figure in the theory part
B_vec = np.linspace(0.0,250.0,1000) # to get the high field limit
B_vec = np.linspace(7,7.5,500) # to get the high field limit
B_vec = np.linspace(7.3,7.42,500) # to get the high field limit
B_vec = np.linspace(6,8,500) # to get the high field limit
# B_vec = np.linspace(7.3,7.362,1000) # below
B_vec = np.linspace(7.362,7.42,1000) # above

# fine grid around 3.629 GHz resonance frequency
# B_vec = np.linspace(0.0822,0.0824,500) # to precisely interpolate the resonance field
# B_vec = np.linspace(0.0815,0.083,500) # to get the deviations on field stepping in rabi oscillations

ge = 28.02495
ge = 1.76085963023e11/2/np.pi*1e-9
w1 = ge*B_vec
gmu = 0.1355

w2 = -gmu*B_vec

# put the exact quantities here
gu =-2.0023318418
uu = -4.49044830e-26
h = 6.62607015e-34

gmu = -gu*uu/h*1e-9

w2 = gmu * B_vec

#J = 4.4956


S = 0.5
I = 0.5


# get all individual operators
Sx_i = qt.jmat(S,"x")
Ix_i = qt.jmat(I,"x")
Sy_i = qt.jmat(S,"y")
Iy_i = qt.jmat(I,"y")
Sz_i = qt.jmat(S,"z")
Iz_i = qt.jmat(I,"z")

Sp_i = qt.jmat(S,"+")
Ip_i = qt.jmat(I,"+")
Sm_i = qt.jmat(S,"-")
Im_i = qt.jmat(I,"-")

Se_i  = qt.qeye(int(2*S+1))
Ie_i  = qt.qeye(int(2*I+1))

a_i = qt.qdiags([1,0],0)
b_i = qt.qdiags([0,1],0)

Sx = qt.tensor(Sx_i, Ie_i)
Sy = qt.tensor(Sy_i, Ie_i)
Sz = qt.tensor(Sz_i, Ie_i)

Sp = qt.tensor(Sp_i, Ie_i)
Sm = qt.tensor(Sm_i, Ie_i)

Ix = qt.tensor(Se_i, Ix_i)
Iy = qt.tensor(Se_i, Iy_i)
Iz = qt.tensor(Se_i, Iz_i)

Ip = qt.tensor(Se_i, Ip_i)
Im = qt.tensor(Se_i, Im_i)

Sz_a = qt.tensor(Sz_i,a_i)
Sz_b = qt.tensor(Sz_i,b_i)
Iz_a = qt.tensor(a_i, Iz_i)
Iz_b = qt.tensor(b_i, Iz_i)

Unity = qt.tensor(Se_i, Ie_i)

psi_mat = np.zeros(Sz.shape)

Fx = Sx.full() + Ix.full()
Fz = Sz.full() + Iz.full()


mixangle = 1/2 * np.arctan2(J,(w1-w2))

Es = []
transs = []
pols = []
# build up the Hamiltonian
for ii_B, B0 in enumerate(B_vec):
    H0_qt = w1[ii_B] * Sz + w2[ii_B] * Iz + J * (Sx*Ix + Sy*Iy + Sz*Iz)
    
    H0 = np.real(H0_qt.full())
    
    psi, E_mat, psi_t = jacobi_math.jacobi_diagonalize(H0)
    
    E = np.diag(E_mat)
    
    # the transformed Fx operator. this will indicate the transition probability
    Fx_t = np.dot(np.dot(psi, Fx),psi_t)
    Ix_t = np.dot(np.dot(psi, Ix.full()),psi_t)
    Sx_t = np.dot(np.dot(psi, Sx.full()),psi_t)
    Iz_t = np.dot(np.dot(psi, Iz.full()),psi_t)
    Fz_t = np.dot(np.dot(psi, Fz),psi_t)
    
    transprob = Fx_t**2
    
    init = (Unity.full() + Iz.full())/4
    init_t = np.dot(np.dot(psi, init),psi_t)
    
    init2 = (Unity.full() + Ix.full())/4
    init2_t = np.dot(np.dot(psi, init2),psi_t)
    
    transprob = np.dot(init2_t,Ix_t)
    
    # excitation Ham in GHz/T
    exc_ham = ge * Sx_t - gmu * Ix_t
    
    # get out the transition energies and intensities
    off_id1, off_id2 = np.triu_indices_from(transprob,k=1)
    
    
    transitions = np.zeros((len(np.triu_indices(int(2*S+1)*int(2*I+1),k=1)[0]),4))
    for jj in range(len(off_id1)):
        transitions[jj,0] = abs(E[off_id1[jj]] - E[off_id2[jj]])
        transitions[jj,1] = np.real(transprob[off_id1[jj], off_id2[jj]])
        transitions[jj,2] = np.real(exc_ham[off_id1[jj], off_id2[jj]])
        transitions[jj,3] = np.real(Iz_t[off_id1[jj], off_id1[jj]]-Iz_t[off_id2[jj], off_id2[jj]])/2
        
    pol = np.zeros((4,1))
    
    for jj in range(len(pol)):
        if jj == 0:
            el = Iz_a
        elif jj == 1:
            el = Iz_b
        elif jj == 2:
            el = Sz_a
        elif jj == 3:
            el = Sz_b
            
        pol[jj] = np.trace(np.dot(Iz_t, el.full()))
    
    Es.append(E)
    transs.append(transitions)
    pols.append(pol)
    
Es = np.array(Es)
transs = np.array(transs)

plt.figure(1)
plt.clf()
plt.plot(B_vec, Es)
plt.ylabel('$E$ [GHz]')

# plt.xlim(B_vec[0],0.5)
# plt.ylim(0,20)
plt.ylabel('Energy [GHz]')
plt.xlabel('Magnetic field $B_0$ [T]')

ii_tr = 5
plt.figure(2)
plt.clf()
plt.subplot(121)
plt.plot(B_vec, transs[:,:, 0])
plt.subplot(122)
plt.plot(B_vec, transs[:,:, 2])

plt.figure(4)
plt.clf()
plt.plot(B_vec, mixangle/np.pi*180)

# plot with variable linewidth
from matplotlib.collections import LineCollection

plt.figure(5)
plt.clf()
a = plt.gca()
# fig,a = plt.subplots()

lwsca = 1.0


# d34_t, d12_t , d13_t, d24_t 
r_ord = [5,0,1,4]
cols=['C0','C1','C2','C3']

r_ord.reverse()
cols.reverse()

# for ii in range(transs.shape[1]):
for ii_c, ii in enumerate(r_ord):
    x = B_vec
    y = transs[:,ii,0].real
    lwidths=transs[:,ii, 2].real*lwsca
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, linewidths=lwidths,color=cols[ii_c])
    
    #
    # a.add_collection(lc)
    # plt.scatter(x,y, s=lwidths**2, color=cols[ii_c])
    plt.plot(x,y*1e3, color=cols[ii_c])
# fig.show()

# plt.plot([0.1,0.2],[10,10])
# plt.xlim(x[0],0.5)
plt.ylim(0,10)
plt.ylim(0,200)
plt.ylabel('$\omega / (2 \pi)$ [MHz]')
plt.xlabel('Magnetic field $B_0$ [T]')

# get partial polarizations at 82.2 mT
t_field = 0.0822
t_idx = np.where(B_vec > t_field)[0][0]

# also get the 12 resonance field at a particular microwave frequency
# requires a finer grid of B around the resonance position to be accurate

t_frq = 0.002
t_frq = 0.004
if t_frq == 0:
    
    t_frq_idx = np.argmin(transs[:,0,0])
else:
    t_frq_idx = np.where(transs[:,0,0] < t_frq)[0][0]

# get 100 points around this (or just half the index in case of a less fine grid)
t_frq_pts = np.min([int(t_frq_idx/4)*2, 100])

t_frqs_idx = np.arange(t_frq_idx - t_frq_pts/2, t_frq_idx + t_frq_pts/2,dtype=int)

srt_ids = np.argsort(transs[t_frqs_idx,0,0])

intp_field = np.interp(t_frq, transs[t_frqs_idx[srt_ids],0,0], B_vec[t_frqs_idx[srt_ids]])

addinfo = 'The grid appears fine'

print('Interpolated resonance field at {} GHz microwave frequency: {:.03f} \n Interpolated from a grid with stepsize of {:.03f} uT'.format(t_frq, intp_field*1e3, np.mean(np.diff(B_vec[t_frqs_idx[srt_ids]]))*1e6))

# also plot the B_0 / w_34 dependence, to check how much one makes an error when linearizing this
plt.figure(6)
plt.clf()
plt.subplot(211)
plt.plot(B_vec*1e3, transs[:,5,0])

# linearize
linpol = np.polyfit(B_vec*1e3, transs[:,5,0], 1);

linfit = np.polyval(linpol, B_vec*1e3)

plt.plot(B_vec*1e3, linfit)

plt.subplot(212)

plt.plot(B_vec*1e3, (transs[:,5,0] - linfit)*1e6)

plt.ylabel('Deviation from linearity [kHz]')

# also plot explicitly the transition moment of the w_34 transition
plt.figure(7)
plt.clf()
plt.plot(B_vec*1e3, transs[:,5,2])
# get the deviation
transmean = np.mean(transs[:,5,2])
devs = [np.min(transs[:,5,2]-transmean), np.max(transs[:,5,2]-transmean)]

#%% additional figure from spidyan sim

import scipy.io
mat = scipy.io.loadmat('Rabi_onres.mat')

sigs = mat['signalc_f'][0][0]
t = mat['t'][0]
scas = mat['d_scalings'][0]
t_rabi = mat['t_rabi'][0][0]

psigs = np.arange(1,sigs.shape[0])

p_idx = t <= 2*t_rabi

p_skip = 1000

for ii, p_id in enumerate(p_idx):
    if (ii%p_skip == 0) and p_id:
        p_idx[ii] = True
    else:
        p_idx[ii] = False
    


plt.figure(6)
plt.clf()
for ii in psigs:
    plt.plot(t[p_idx]/t_rabi, sigs[ii,p_idx]*scas[ii-1])

plt.plot(t[p_idx]/t_rabi, np.sum(sigs[1:,p_idx]*scas[:,None],axis=0))

plt.xlabel('$t / T_\mathrm{Rabi}$')
plt.ylabel('$<I_z>$')



#%% yet another extrafig from spidyan sims with offresonance

# save('Rabi_offres','fitpars','W_vec', 't_rabi','Iz_norm','d_scalings');
mat = scipy.io.loadmat('Rabi_offres.mat')


fitpars = mat['fitpars']
fitpars_full = mat['fitpars_full']
W_vec = mat['W_vec'][0]
t_rabi = mat['t_rabi'][0]
Iz_norm = mat['Iz_norm'][0]
d_scalings = mat['d_scalings'][0]

# the partial polarizations: square matrix elements (with appropriate normalization)
partpol = d_scalings**2 * 2

nu1 = 1/t_rabi


weff_v = np.sqrt((nu1)**2 + W_vec**2);
amps = partpol[0] * (1.0 - abs(W_vec/weff_v)**2);

plt.figure(7)
plt.clf()
plt.plot(W_vec*1e3, fitpars[0,:])
plt.plot(W_vec*1e3, fitpars_full[0,:])
plt.plot(W_vec*1e3, amps)
plt.xlabel('$\Omega / (2\pi)$ [MHz]')
plt.ylabel('$\Delta$$A$')

plt.figure(8)
plt.clf()
plt.plot(W_vec*1e3, fitpars[1,:]*1e3)
plt.plot(W_vec*1e3, fitpars_full[1,:]*1e3)
plt.plot(W_vec*1e3, weff_v*1e3)
plt.xlabel('$\Omega / (2\pi)$ [MHz]')
plt.ylabel('$\omega_\mathrm{eff} / (2\pi)$ [MHz]')


plt.figure(9)
plt.clf()
# plt.plot(W_vec*1e3, fitpars[2,:])
plt.plot(W_vec*1e3, fitpars_full[2,:])
plt.plot(W_vec*1e3, np.sum(partpol) - amps)
plt.xlabel('$\Omega / (2\pi)$ [MHz]')
plt.ylabel('$A_\mathrm{bg}$')


plt.figure(10)
plt.clf()
plt.plot(W_vec/nu1, fitpars[2,:])
plt.plot(W_vec/nu1, fitpars_full[2,:])
plt.plot(W_vec/nu1, np.sum(partpol) - amps)
plt.xlabel('$\Omega / (2\pi)$ [MHz]')
plt.ylabel('$A_0$')





#%% some numerical checks to assure that the analytical rules are fulfilled

# diagonlization operator
init = Sz-Iz
init_m = init.full()


fin = 2*Ix*Sx+2*Sy*Iy
fin_m = fin.full()

rot_commutator = np.matmul(init_m,fin_m) - np.matmul(fin_m, init_m)

rot = 2*Sx*Iy - 2*Sy*Ix
rot_m = rot.full()

# first check for the last calculated point
rot_ang = mixangle[-1]

unitary = linalg.expm(-1j*rot_ang*rot_m)
unitary_t = linalg.expm(1j*rot_ang*rot_m)

H_diag = np.matmul(np.matmul(unitary,H0),unitary_t)



Fx_t2 = np.matmul(np.matmul(unitary, Fx),unitary_t)

# eigvec = np.matmul(np.matmul(unitary,np.eye(4)),unitary_t)


print('Numerical eigenvalues')
print(E_mat.real)

print('Analytical eigenvalues')
print(H_diag.real)



print('Numerical eigenvectors')
print(psi.real)

print('Analytical unitary')
print(unitary.real)




print('Numerical Fx')
print(Fx_t.real)

print('Analytical Fx')
print(Fx_t2.real)



