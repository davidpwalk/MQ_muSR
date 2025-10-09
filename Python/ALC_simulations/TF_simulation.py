import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import qutip as qt
import xarray as xr
from scipy.fft import fft, fftfreq, fftshift

from Python.jacobi_math import jacobi_diagonalize

# Set custom plotting style
from Python.plot_settings import set_demonlab_style
set_demonlab_style()
pio.templates.default = "DemonLab"
#%% Parameters and calculations of spin operators and tensors
# Magnetic field range (in Tesla)
magnetic_fields = np.linspace(0, 5, 100)
# magnetic_fields = [0]

# Zeeman (gamma / GHz/T)
ge = 28.02495
gmu = -0.1355

# Coupling constants (in GHz) and rotation angle (in degrees)
# A_iso = 0.5148
A_iso = 0.0
D_parallel = 0.002
# D_parallel = 0
D_perp = -D_parallel/2

# Rotation angles (in degrees)
thetas = np.radians(np.linspace(0, 90, 200, dtype=np.float64))
# thetas = np.radians([0, 45, 90, 180])
# thetas = [0]

# Define the spin operators for S=1/2 and I=1/2
S = 0.5
I = 0.5

Sx = qt.tensor(qt.jmat(S, 'x'), qt.qeye(int(2*I + 1)))
Sy = qt.tensor(qt.jmat(S, 'y'), qt.qeye(int(2*I + 1)))
Sz = qt.tensor(qt.jmat(S, 'z'), qt.qeye(int(2*I + 1)))
Ix = qt.tensor(qt.qeye(int(2*S + 1)), qt.jmat(I, 'x'))
Iy = qt.tensor(qt.qeye(int(2*S + 1)), qt.jmat(I, 'y'))
Iz = qt.tensor(qt.qeye(int(2*S + 1)), qt.jmat(I, 'z'))

S_ops = [Sx, Sy, Sz]
I_ops = [Ix, Iy, Iz]

T_principal = np.diag([D_perp, D_perp, D_parallel])


# Rotation matrix
def Ry(theta):
    return np.array([[np.cos(theta),  0, np.sin(theta)],
                     [0,              1,             0],
                     [-np.sin(theta), 0, np.cos(theta)]])

# Generate rotated tensors
T_labs = []
for theta in thetas:
    R = Ry(theta)
    T_rot = R @ T_principal @ R.T
    T_labs.append(T_rot)

#%% Simulation
# Settings
O = Ix+Sx  # observable
threshold = 1e-4  # amplitude threshold for transitions
transition_type_filter = None  # "SQE", "SQMu",  "ZQ", "DQ" or None for all


time = np.linspace(0, 8000, 16000)
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
    elif pair == {1, 2} or pair == {2, 4}:
        return "SQE"
    elif pair == {1, 3} or pair == {3, 4}:
        return "SQMu"


import numpy as np

def merge_transitions(freqs, amps, types, tol=1e-8, zero_tol=1e-12):
    """
    Merge transitions whose frequencies are pairwise within `tol` (single-linkage / connected components).
    - freqs, amps, types: sequences (will be converted to numpy arrays)
    - tol: absolute frequency tolerance for linking two transitions
    - zero_tol: if set, snap |freq| <= zero_tol to 0 before clustering (useful for near-zero symmetric pairs)
    Returns: (merged_freqs, merged_amps, merged_types) as Python lists (sorted by frequency).
    """
    freqs = np.asarray(freqs, dtype=float)
    amps = np.asarray(amps, dtype=float)
    types = np.asarray(types, dtype=object)

    if freqs.size == 0:
        return [], [], []

    # optional: snap tiny frequencies to zero
    if zero_tol is not None:
        small = np.abs(freqs) <= zero_tol
        freqs[small] = 0.0

    # adjacency: True where pairwise difference <= tol
    adj = np.abs(freqs[:, None] - freqs[None, :]) <= tol

    # find connected components (DFS)
    n = len(freqs)
    visited = np.zeros(n, dtype=bool)
    groups = []
    for i in range(n):
        if visited[i]:
            continue
        stack = [i]
        comp = []
        visited[i] = True
        while stack:
            v = stack.pop()
            comp.append(v)
            neighbors = np.where(adj[v] & ~visited)[0]
            for w in neighbors:
                visited[w] = True
                stack.append(w)
        groups.append(comp)

    # merge each component: amplitude sum, frequency = amplitude-weighted mean (fallback to simple mean)
    merged_freqs = []
    merged_amps = []
    merged_types = []
    for comp in groups:
        idx = np.array(comp)
        total_amp = amps[idx].sum()
        if total_amp > 0:
            freq_rep = (freqs[idx] * amps[idx]).sum() / total_amp
        else:
            freq_rep = freqs[idx].mean()
        merged_freqs.append(freq_rep)
        merged_amps.append(total_amp)
        merged_types.append(",".join(sorted(set(types[idx].tolist()))))

    # sort by frequency before returning
    order = np.argsort(merged_freqs)
    merged_freqs = np.array(merged_freqs)[order].tolist()
    merged_amps  = np.array(merged_amps)[order].tolist()
    merged_types = np.array(merged_types, dtype=object)[order].tolist()

    return merged_freqs, merged_amps, merged_types


for ii, T_lab in enumerate(T_labs):
    for jj, magnetic_field in enumerate(magnetic_fields):
        nu_muon = gmu * magnetic_field
        nu_electron = ge * magnetic_field

        H_zeeman = nu_electron * Sz + nu_muon * Iz
        H_iso = A_iso * (Sx * Ix + Sy * Iy + Sz * Iz)
        H_dip = qt.Qobj(np.zeros(Sx.shape, dtype=complex), dims=Sx.dims)
        for a in range(len(S_ops)):
            for b in range(len(S_ops)):
                scalar = float(T_lab[a, b])
                H_dip += scalar * (S_ops[a] * I_ops[b])

        # H_dip_approx = D_parallel * Sz * Iz + D_perp * (Sx * Ix + Sy * Iy)

        H_tot = H_zeeman + H_iso + H_dip

        # Diagonalize the Hamiltonian using the Jacobi method (H_tot needs to be converted from Qobj to ndarray)
        psi, H_tot_eigen, psi_t = jacobi_diagonalize(np.real(H_tot.full()))

        energies = H_tot_eigen.diagonal()
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
                # print(f'ttype: {ttype}, kk: {kk+1}, ll: {ll+1}, amp: {amp}, O: \n{O_eigen}') if ttype == 'DQ' or ttype == 'ZQ' else None
                print(f'theta: {np.degrees(thetas[ii])}, energies: {energies}, \n T_lab: {T_lab} \nH_tot_eigen: \n{H_tot_eigen}')
                if amp > threshold:
                    freqs.append(energies[ll] - energies[kk])
                    amps.append(amp)
                    ttype = classify_transition(ll + 1, kk + 1)
                    types.append(ttype)

        # store frequencies/amplitudes (pad with NaN if fewer than max_transitions)
        freqs, amps, types = merge_transitions(freqs, amps, types)
        print(f'freqs: {freqs}, amps: {amps}, types: {types}')
        n_trans = len(freqs)
        frequencies_arr[ii, jj, :n_trans] = freqs
        amplitudes_arr[ii, jj, :n_trans] = amps
        ttypes_arr[ii, jj, :n_trans] = types

        # build signal
        signal = np.zeros_like(time, dtype=float)
        for amp, freq in zip(amps, freqs):
            signal += amp * np.cos(2 * np.pi * abs(freq) * time)

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
        "transition_type": ["SQE", "SQMu", "ZQ", "DQ"],
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
    else:
        freqs = results['frequencies'].sel(theta=theta, B=B, method='nearest').values
        amps = results['amplitudes'].sel(theta=theta, B=B, method='nearest').values

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
        freqs = abs(results['frequencies'].sel(theta=theta, B=B, method='nearest').where(mask, drop=True).values)
        amps = results['amplitudes'].sel(theta=theta, B=B, method='nearest').where(mask, drop=True).values
        ttypes = results['transition_types'].sel(theta=theta, B=B, method='nearest').where(mask, drop=True).values
    else:
        freqs = abs(results['frequencies'].sel(theta=theta, B=B, method='nearest').values)
        amps = results['amplitudes'].sel(theta=theta, B=B, method='nearest').values
        ttypes = results['transition_types'].sel(theta=theta, B=B, method='nearest').values

    # remove NaN padding
    mask = np.isfinite(freqs) & np.isfinite(amps)
    freqs = freqs[mask]
    amps = amps[mask]
    idx = 0

    fig = go.Figure()
    for freq, amp, ttype in zip(freqs, amps, ttypes):
        color = px.colors.qualitative.G10[idx]
        fig.add_trace(
            go.Scatter(
                x=[freq, freq],
                y=[0, amp],
                mode="lines",
                name=f"{freq:.4f} GHz; {amp:.4f}, {ttype}",
                line=dict(color=color),
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[freq], y=[amp],
                mode="markers",
                marker=dict(color="rgba(0,0,0,0)"),
                showlegend=False,
            )
        )
        idx = idx + 1

    fig.add_hline(y=0, line=dict(color="black"), opacity=1)

    fig.update_layout(
        xaxis_title="Frequency / GHz",
        yaxis_title="Amplitude",
        template="DemonLab"
    )
    return fig

#%% Plotting single spectra examples
fig = time_signal(0, 0)
fig.update_layout(xaxis_range=[0, 30])
# fig.show()

fig = stick_spectrum(0, 0)
fig.update_layout(title='TF muSR spectrum at θ = 0°, B = 0 T, O = Sx + Ix')
fig.show()
#%% Plot angular dependence of spectra at fixed B
B = 1  # Tesla
thetas = np.linspace(0, 90, 6)
colors = px.colors.qualitative.G10[:len(thetas)]

# Create a common figure
fig = go.Figure()

for theta, color in zip(thetas, colors):
    subfig = stick_spectrum(theta, B=B)
    for trace in subfig.data:
        # Only modify style of the line traces
        if trace.mode == "lines":
            trace.line.color = color
            trace.showlegend = False
        else:
            trace.marker.color = color
            trace.name = f"{theta}°"
        fig.add_trace(trace)

fig.add_hline(y=0, line=dict(color="black"), opacity=1)

fig.update_layout(
    xaxis_title="Frequency / GHz",
    yaxis_title="Amplitude",
    title=f'Angular dependence of muSR spectrum at B = {B} T',
    template="DemonLab",
    legend=dict(y=0.5, x=1.02, xanchor="left", yanchor="middle"),
)
fig.show(renderer='browser')
# fig.write_html('../../Figures/ALC_simulations/TF_angular_dependence.html')

#%% Plot B dependence of spectra at fixed angle
magnetic_fields = np.linspace(0, 1, 6)
colors = px.colors.qualitative.G10[:len(magnetic_fields)]

# Create a common figure
fig = go.Figure()

for magnetic_field, color in zip(magnetic_fields, colors):
    subfig = stick_spectrum(0, B=magnetic_field)
    for trace in subfig.data:
        # Only modify style of the line traces
        if trace.mode == "lines":
            trace.line.color = color
            trace.showlegend = False
        else:
            trace.marker.color = color
            trace.name = f"{magnetic_field:.1f} T"
        fig.add_trace(trace)

fig.add_hline(y=0, line=dict(color="black"), opacity=1)

fig.update_layout(
    xaxis_title="Frequency / GHz",
    yaxis_title="Amplitude",
    title='B dependence of muSR spectrum at θ = 0',
    template="DemonLab",
    legend=dict(y=0.5, x=1.02, xanchor="left", yanchor="middle"),
)
fig.show(renderer='browser')
# fig.write_html('../../Figures/ALC_simulations/TF_B_dependence.html')
#%% Sum time-domain signals over angles
# TODO: check impact of weighting and normalization
# powder_signals = signals.sum(dim='theta')

weights = np.sin(np.radians(signals.theta))
powder_signals = (signals * weights).sum(dim='theta') / weights.sum()

fig = px.line(x=powder_signals.sel(B=1, method='nearest').time, y=powder_signals.sel(B=1, method='nearest'), labels={'x': 'Time / ns',
                                                                                                                     'y': 'Signal'})
fig.show()

def apodize(signal):
    n = len(signal)
    dt = time[1] - time[0]
    window = np.hanning(2*n)[n:]
    apodized_signal = signal * window

    n_fft = 1 << (n - 1).bit_length()  # next power of 2
    padded = np.zeros(n_fft)
    padded[:n] = apodized_signal
    time_padded = np.arange(n_fft) * dt
    return padded, time_padded

def ft(signal):
    signal, time = apodize(signal)
    fig = px.line(x=time, y=signal, labels={'x': 'Time / ns', 'y': 'Signal'})
    fig.show()
    dt = time[1] - time[0]
    power_spectrum = (np.abs(fft(signal)))**2
    power_spectrum = power_spectrum/max(power_spectrum)
    freq = fftshift(fftfreq(len(signal), d=dt))
    return power_spectrum, freq

def plot_powder_spectrum(signal):
    spectrum, freq = ft(signal)
    fig = px.line(x=freq, y=spectrum, labels={'x': 'Frequency / GHz', 'y': 'Intensity'})
    # fig = px.line(x=np.arange(spectrum.shape[0]), y=spectrum, labels={'x': 'Frequency / GHz', 'y': 'Intensity'})
    return fig

#%% Plot powder spectra
B = 0.05  # Tesla
fig = plot_powder_spectrum(powder_signals.sel(B=B, method='nearest').values)
fig.update_yaxes(automargin=True)
fig.update_layout(title=f'Powder muSR spectrum at B = {B} T, O = Ix + Sx')
fig.show()
fig.write_html(f'../../Figures/ALC_simulations/TF_powder_spectrum_thin_pake_pattern_{B}T_Ix_Sx.html')