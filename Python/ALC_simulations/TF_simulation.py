import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import qutip as qt
import xarray as xr
from scipy.fft import fft, fftfreq

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
# A_iso = 0.5148
A_iso = 0.0
D_parallel = 0.002
D_parallel = 0.3
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
# Settings
O = Ix  # observable
threshold = 1e-6  # amplitude threshold for transitions
transition_type_filter = None  # "SQ", "ZQ", "DQ" or None for all


time = np.linspace(0, 500, 2000)
n_theta = len(thetas)
n_B = len(magnetic_fields)
n_time = len(time)

# Transition count is small: max is 4 levels → 6 transitions
max_transitions = 6

energies_arr = np.zeros((n_theta, n_B, 4))  # 4 levels for S=I=1/2
frequencies_arr = np.full((n_theta, n_B, max_transitions), np.nan)
amplitudes_arr = np.full((n_theta, n_B, max_transitions), np.nan)
ttypes_arr = np.full((n_theta, n_B, max_transitions), None, dtype=object)

signals_arr = np.zeros((n_theta, n_B, n_time))

def classify_transition(i, j):
    pair = {i, j}
    if pair == {1, 4}:
        return "DQ"
    elif pair == {2, 3}:
        return "ZQ"
    else:
        return "SQ"

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

        O_dense = O.full()
        O_eigen = psi.conj().T @ O_dense @ psi

        freqs = []
        amps = []
        types = []
        for kk in range(len(energies)):
            for ll in range(kk + 1, len(energies)):
                if transition_type_filter:
                    ttype = classify_transition(ll + 1, kk + 1)
                    if ttype != transition_type_filter:
                        continue

                amp = abs(O_eigen[kk, ll]) ** 2
                ttype = classify_transition(ll + 1, kk + 1)
                print(f'ttype: {ttype}, kk: {kk + 1}, ll: {ll + 1}, amp: {amp}, O: {O_dense}') if ttype == "DQ" or ttype == "ZQ" else None
                if amp > threshold:
                    freqs.append(energies[ll] - energies[kk])
                    amps.append(amp)

                    ttype = classify_transition(ll + 1, kk + 1)
                    types.append(ttype)

        # store frequencies/amplitudes (pad with NaN if fewer than max_transitions)
        n_trans = len(freqs)
        frequencies_arr[ii, jj, :n_trans] = freqs
        amplitudes_arr[ii, jj, :n_trans] = amps
        ttypes_arr[ii, jj, :n_trans] = types

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
        "transition_types": (["theta", "B", "transition"], ttypes_arr),
    },
    coords={
        "theta": np.degrees(thetas),
        "B": magnetic_fields,
        "time": time,
        "level": np.arange(4),
        "transition": np.arange(max_transitions),
        "transition_type": ["SQ", "ZQ", "DQ"],
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

#%% Plotting functions for single spectra
def time_signal(theta, B, transition_type=None):
    if transition_type:
        mask = results['transition_types'].sel(theta=theta, B=B, method='nearest') == transition_type
        freqs = results['frequencies'].sel(theta=theta, B=B, method='nearest').where(mask, drop=True).values
        amps = results['amplitudes'].sel(theta=theta, B=B, method='nearest').where(mask, drop=True).values

        # remove NaN padding
        mask = np.isfinite(freqs) & np.isfinite(amps)
        freqs = freqs[mask]
        amps = amps[mask]

        signal = np.zeros_like(time, dtype=float)
        for amp, freq in zip(amps, freqs):
            signal += amp * np.cos(2 * np.pi * freq * time)

        fig = px.line(x=time, y=signal, labels={'x': 'Time / ns', 'y': 'Signal'})
        return fig

def stick_spectrum(theta, B, transition_type=None):
    if transition_type:
        mask = results['transition_types'].sel(theta=theta, B=B, method='nearest') == transition_type
        freqs = results['frequencies'].sel(theta=theta, B=B, method='nearest').where(mask, drop=True).values
        amps = results['amplitudes'].sel(theta=theta, B=B, method='nearest').where(mask, drop=True).values
    else:
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
                name=f"{f:.3f} GHz; {a:.3f}, {transition_type}",
            )
        )
    if len(freqs) > 1:
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

#%% Plotting single spectra examples
# fig = time_signal(4, 6)
# fig.show()
#
fig = stick_spectrum(4, 6, transition_type="SQ")
fig.show()

#%% Sum time-domain signals over angles
# TODO: check impact of weighting and normalization
# powder_signals = signals.sum(dim='theta')

weights = np.sin(np.radians(signals.theta))
powder_signals = (signals * weights).sum(dim='theta') / weights.sum()

fig = px.line(x=powder_signals.sel(B=2).time, y=powder_signals.sel(B=6))
fig.show()

def apodize(signal):
    n = len(signal)
    window = np.hanning(2*n)[n:]
    apodized_signal = signal * window

    n_fft = 1 << (n - 1).bit_length()  # next power of 2
    padded = np.zeros(n_fft)
    padded[:n] = apodized_signal
    return padded

def ft(signal):
    signal = apodize(signal)
    power_spectrum = (abs(fft(signal))**2)[:len(signal)//2]
    freq = fftfreq(len(signal), d=1/len(signal))[:len(signal)//2]
    return power_spectrum, freq

def plot_powder_spectrum(signal):
    spectrum, freq = ft(signal)
    fig = px.line(x=freq, y=spectrum, labels={'x': 'Frequency / GHz', 'y': 'Intensity'})
    fig.show()
    return fig

#%% Plot powder spectra

powder_spectrum, freq = ft(powder_signals.sel(B=2))
fig = px.line(x=freq, y=powder_spectrum, labels={'x': 'Frequency / GHz', 'y': 'Intensity'})
fig.show()