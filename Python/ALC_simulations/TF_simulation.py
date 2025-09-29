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
# thetas = np.radians([45])

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
n_theta = len(thetas)
n_B = len(magnetic_fields)
n_time = len(time)

# Transition count is small: max is 4 levels → 6 transitions
max_transitions = 6

energies_arr = np.zeros((n_theta, n_B, 4))  # 4 levels for S=I=1/2
frequencies_arr = np.full((n_theta, n_B, max_transitions), np.nan)
amplitudes_arr = np.full((n_theta, n_B, max_transitions), np.nan)

signals_arr = np.zeros((n_theta, n_B, n_time))

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
        energies_arr[ii, jj, :] = energies

        Ix_dense = Ix.full()
        Ix_eigen = psi.conj().T @ Ix_dense @ psi

        freqs = []
        amps = []
        for kk in range(len(energies)):
            for ll in range(kk + 1, len(energies)):
                freqs.append(energies[ll] - energies[kk])
                amps.append(abs(Ix_eigen[kk, ll]) ** 2)

        # store frequencies/amplitudes (pad with NaN if fewer than max_transitions)
        n_trans = len(freqs)
        frequencies_arr[ii, jj, :n_trans] = freqs
        amplitudes_arr[ii, jj, :n_trans] = amps

        # build signal
        signal = np.zeros_like(time, dtype=float)
        for amp, freq in zip(amps, freqs):
            signal += amp * np.cos(2 * np.pi * freq * time)

        signals_arr[ii, jj, :] = signal

# Wrap into Dataset
results = xr.Dataset(
    {
        "signal": (["theta", "B", "time"], signals_arr),
        "energies": (["theta", "B", "level"], energies_arr),
        "frequencies": (["theta", "B", "transition"], frequencies_arr),
        "amplitudes": (["theta", "B", "transition"], amplitudes_arr),
    },
    coords={
        "theta": np.degrees(thetas),
        "B": magnetic_fields,
        "time": time,
        "level": np.arange(4),
        "transition": np.arange(max_transitions),
    }
)

signals = xr.DataArray(
    signals_arr,
    dims=["theta", "B", "time"],
    coords={
        "theta": np.degrees(thetas),
        "B": magnetic_fields,
        "time": time,
    },
    name="signal")

#%% Plotting Functions
def time_signal(theta, B):
    signal = signals.sel(theta=theta, B=B, method='nearest')
    fig = px.line(x=signal.time, y=signal, labels={'x': 'Time / ns', 'y': 'Signal'})
    return fig

def stick_spectrum(theta, B):
    freqs = results['frequencies'].sel(theta=theta, B=B, method='nearest').values
    amps = results['amplitudes'].sel(theta=theta, B=B, method='nearest').values

    # remove NaN padding
    mask = np.isfinite(freqs) & np.isfinite(amps)
    freqs = freqs[mask]
    amps = amps[mask]

    fig = go.Figure()
    for f, a in zip(freqs, amps):
        fig.add_trace(
            go.Scatter(
                x=[f, f],
                y=[0, a],
                mode="lines",
                line=dict(color="black"),
                showlegend=False
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[f], y=[a],
                mode="markers",
                marker=dict(color="rgba(0,0,0,0)"),
                name=f"{f:.3f} GHz; {a:.3f}"
            )
        )

    fig.add_trace(
        go.Scatter(
            x=[freqs.min() - 1, freqs.max() + 1],
            y=[0, 0],
            mode="lines",
            line=dict(color="black"),
            showlegend=False
        )
    )

    fig.update_layout(
        xaxis_title="Frequency / GHz",
        yaxis_title="Amplitude",
        template="DemonLab"
    )
    return fig

#%% Plotting
fig = time_signal(4, 0)
fig.show()

fig = stick_spectrum(4, 5)
fig.show()