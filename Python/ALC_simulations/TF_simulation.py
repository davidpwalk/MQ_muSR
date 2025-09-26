import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import qutip as qt
import xarray as xr

from Python.jacobi_math import jacobi_diagonalize

# Set custom plotting style
from Python.plot_settings import set_demonlab_style
set_demonlab_style()
pio.templates.default = "DemonLab"
#%% Parameters and calculations of spin operators and tensors
# Magnetic field range (in Tesla)
magnetic_fields = np.linspace(0, 12, 13)

# Zeeman (gamma / GHz/T)
ge = 28.02495
gmu = -0.1355

# Coupling constants (in GHz) and rotation angle (in degrees)
A_iso = 0.5148
D_parallel = 0.002
D_perp = -D_parallel/2

# Rotation angles (in degrees)
thetas = np.radians(np.linspace(0, 90, 400, dtype=np.float64))
thetas = np.radians([45])

# Define the spin operators for S=1/2 and I=1/2
S = 0.5
I = 0.5

Sx = qt.tensor(qt.jmat(S, 'x'), qt.qeye(int(2*I + 1)))
Sy = qt.tensor(qt.jmat(S, 'y'), qt.qeye(int(2*I + 1)))
Sz = qt.tensor(qt.jmat(S, 'z'), qt.qeye(int(2*I + 1)))
Ix = qt.tensor(qt.qeye(int(2*S + 1)), qt.jmat(I, 'x'))
Iy = qt.tensor(qt.qeye(int(2*S + 1)), qt.jmat(I, 'y'))
Iz = qt.tensor(qt.qeye(int(2*S + 1)), qt.jmat(I, 'z'))

T_principal = np.diag([D_perp, D_perp, D_parallel])


# Rotation matrix
# TODO: check if this is the right one to remain (is Ry or Rz?)
def Rz(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta),  np.cos(theta), 0],
                     [0,            0,           1]])

# Generate rotated tensors (T_labs are of type Qobj)
T_labs = []
for theta in thetas:
    R = Rz(theta)
    T_rot = R @ T_principal @ R.T
    T_labs.append(qt.Qobj(T_rot))

#%% Simulation
time = np.linspace(0, 500, 2000)
results = []
signals = np.zeros((len(T_labs), len(magnetic_fields), len(time)))

for ii, T_lab in enumerate(T_labs):
    for jj, magnetic_field in enumerate(magnetic_fields):
        nu_muon = gmu * magnetic_field
        nu_electron = ge * magnetic_field

        H_zeeman = nu_electron * Sz + nu_muon * Iz
        H_iso = A_iso * (Sx * Ix + Sy * Iy + Sz * Iz)
        H_dip = (Sx * (T_lab[0, 0] * Ix + Sx * T_lab[0, 1] * Iy + Sx * T_lab[0, 2] * Iz) +
                 Sy * (T_lab[1, 0] * Ix + Sy * T_lab[1, 1] * Iy + Sy * T_lab[1, 2] * Iz) +
                 Sz * (T_lab[2, 0] * Ix + Sz * T_lab[2, 1] * Iy + Sz * T_lab[2, 2] * Iz))

        H_tot = H_zeeman + H_iso + H_dip

        # Diagonalize the Hamiltonian using the Jacobi method (H_tot needs to be converted from Qobj to ndarray)
        psi, energy_matrix, psi_t = jacobi_diagonalize(np.real(H_tot.full()))

        energies = energy_matrix.diagonal()

        # transform observable into eigenbasis
        Ix_dense = Ix.full()
        Ix_eigen = psi.conj().T @ Ix_dense @ psi

        frequencies = []
        amplitudes = []

        for kk in range(len(energies)):
            for ll in range(kk+1, len(energies)):
                freq = energies[ll] - energies[kk]
                amp = abs(Ix_eigen[kk, ll])**2
                frequencies.append(freq)
                amplitudes.append(amp)

        signal = np.zeros_like(time, dtype=float)
        for amp, freq in zip(amplitudes, frequencies):
            signal += amp * np.cos(2 * np.pi * freq * time)

        results.append({
            'theta_idx': ii,
            'theta': np.degrees(thetas[ii]),
            'magnetic_field': magnetic_field,
            'signal': signal,
            'frequencies': frequencies,
            'amplitudes': amplitudes,
            'H_tot': H_tot,
            'energies': energies
        })

        signals[ii, jj, :] = signal

# Convert results to xarray DataArray for easier handling
signals = xr.DataArray(
    signals,
    dims=["theta", "B", "time"],
    coords={
        "theta": np.degrees(thetas),
        "B": magnetic_fields,
        "time": time,
    },
    name="signal"
)

#%% Plotting
# Time-domain signal
fig = px.line(x=time, y=signals.sel(theta=45, B=4))
fig.show()