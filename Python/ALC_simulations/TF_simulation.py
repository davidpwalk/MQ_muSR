import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import qutip as qt
import xarray as xr
import os
import warnings

from Python.diagonalizer import diagonalize_ordered

# Import utility functions
from Python.ALC_simulations.utils import *

# Set custom plotting style
from Python.plot_settings import set_demonlab_style
set_demonlab_style()
pio.templates.default = "DemonLab"

#%% Parameters and calculations of spin operators and tensors
# Filename and description for saving results (set filename to None to skip saving)
# filename = 'Python/ALC_simulations/Data/TF_NMR_isotropic_cc.nc'
filename = None
desc = 'TF NMR simulation for S=1/2, I=1/2 system with isotropic hyperfine coupling of A_iso=2 MHz, O=Ix'

# TODO: maybe implement option to generate all signals (for each B and theta), this consumes ALOT of memory
gen_all_signals = False

# Magnetic field range (in Tesla)
# magnetic_fields = np.linspace(0.1, 0, 400)
# magnetic_fields = [0, 0.01, 5]
magnetic_fields = [10]

# Zeeman (gamma / GHz/T)
ge = -28.02495
gmu = 0.1355

# Coupling constants (in GHz) and rotation angle (in degrees)
# A_iso = 0.5148
A_iso = 0
# D_parallel = 0.002
D_parallel = 0.005
D_perp = -D_parallel/2

# Rotation angles (in degrees)
thetas_deg = np.linspace(0, 90, 1600, dtype=np.float64)
# thetas_deg = [0, 45, 90]
# thetas_deg = [45]
thetas = np.radians(thetas_deg)

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
O = Sx  # observable
O_string = 'Sx'
threshold = 1e-4  # amplitude threshold for transitions

time = np.linspace(0, 8000, 32000)
n_theta = len(thetas)
n_B = len(magnetic_fields)
n_time = len(time)

# Transition count is small: max is 4 levels → 6 transitions
max_transitions = 6

energies_arr = np.zeros((n_theta, n_B, 4))  # 4 levels for S=I=1/2
frequencies_arr = np.full((n_theta, n_B, max_transitions), np.nan)
amplitudes_arr = np.full((n_theta, n_B, max_transitions), np.nan)
ttypes_arr = np.full((n_theta, n_B, max_transitions), None, dtype=object)

transition_types = ["SQMu_a", "SQMu_b", "SQE_a", "SQE_b", "ZQ", "DQ"]

for ii, T_lab in enumerate(T_labs):
    prev_vecs = None  # for tracking eigenvector continuity (reset for each theta)

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

        if prev_vecs is None:
            prev_vecs = np.eye(H_tot.shape[0], dtype=complex)

        # Diagonalize the Hamiltonian using the Jacobi method (H_tot needs to be converted from Qobj to ndarray)
        energies, psi = diagonalize_ordered(np.real(H_tot.full()), prev_vecs=prev_vecs)

        energies_arr[ii, jj, :] = energies
        prev_vecs = psi

        O_dense = O.full()
        O_eigen = psi.conj().T @ O_dense @ psi

        freqs = []
        amps = []
        ttypes = []
        # print(f'O: \n{O_eigen}')
        for kk in range(len(energies)):
            for ll in range(kk + 1, len(energies)):
                amp = abs(O_eigen[kk, ll]) ** 2
                ttype = classify_transition(kk + 1, ll + 1)
                # print(f'ttype: {ttype}, kk: {kk + 1}, ll: {ll + 1}, amp: {amp}')
                # print(f'ttype: {ttype}, kk: {kk+1}, ll: {ll+1}, amp: {amp}, O: \n{O_eigen}') if ttype == 'DQ' or ttype == 'ZQ' else None
                # print(f'theta: {np.degrees(thetas[ii])}, energies: {energies}, \n T_lab: {T_lab}')
                # print(ttype, kk+1, ll+1)
                ttype = classify_transition(ll + 1, kk + 1)
                ttypes.append(ttype)
                if amp > threshold:
                    freqs.append(energies[ll] - energies[kk])
                    amps.append(amp)
                else:
                    freqs.append(np.nan)
                    amps.append(np.nan)

        # store frequencies/amplitudes (pad with NaN if fewer than max_transitions)
        # freqs, amps, ttypes = merge_transitions(freqs, amps, ttypes)
        # print(f'freqs: {freqs}, amps: {amps}, types: {ttypes}')
        n_trans = len(freqs)
        frequencies_arr[ii, jj, :n_trans] = freqs
        amplitudes_arr[ii, jj, :n_trans] = amps
        ttypes_arr[ii, jj, :n_trans] = ttypes

results = xr.Dataset(
    {
        "energies": (["theta", "B", "level"], energies_arr),
        "frequencies": (["theta", "B", "transition"], frequencies_arr),
        "amplitudes": (["theta", "B", "transition"], amplitudes_arr),
        "transition_types": (["theta", "B", "transition"], ttypes_arr),
    },
    coords={
        "transition_type": transition_types,
        "theta": thetas_deg,
        "B": magnetic_fields,
        "time": time,
        "level": np.arange(4),
        "transition": np.arange(max_transitions),
    }
)

# Attach simulation parameters as metadata
results.attrs.update({
    "description": desc,
    "A_iso_GHz": A_iso,
    "D_parallel_GHz": D_parallel,
    "D_perp_GHz": D_perp,
    "ge_GHz_per_T": ge,
    "gmu_GHz_per_T": gmu,
    "observable": O_string,
})

if filename:
    if os.path.exists(filename):
        warnings.warn(f"File '{filename}' already exists. Results not overwritten.", UserWarning)
    else:
        results.to_netcdf(filename, engine='netcdf4')
        print(f"Results successfully saved to '{filename}'.")
else:
    print("No filename provided. Results not saved to file.")

#%% Plotting single spectra
B = magnetic_fields[0]  # Tesla
theta = np.degrees(thetas[0])  # degrees
# fig = time_signal(results, theta, B)
# fig.update_layout(xaxis_range=[0, 30])
# fig.show()

fig = stick_spectrum(results, theta, B, transition_type=None)
fig.update_layout(
    title=f'TF muSR @ θ = {theta:.0f}°, B = {B:.0f} T, O = {O_string}',
    margin=dict(t=75),
    showlegend=True,
)
fig.show()
print(os.getcwd())
fig.write_html(f'Figures/TF_simulations/TF_ZF_O{O_string}_A{A_iso}_T{D_parallel}.html')

#%% Calculate and plot powder signals
B = 10  # Tesla
# B = magnetic_fields[99]
print(f'B = {B} T')
transition_filter = None
# transition_filter = ['ZQ', 'DQ']
# transition_filter = ['SQMu_a', 'SQMu_b', 'SQE_a', 'SQE_b']
# transition_filter = ['SQE_a']

powder_signals_df = generate_powder_signals(results, time, magnetic_field=B, transition_filter=transition_filter, make_plot=True)

powder_signal = powder_signals_df['Total'].values
time = powder_signals_df.index.values

# px.line(x=time, y=powder_signal).show()

# #%% Plot powder spectra
fig, freq = plot_powder_spectrum(powder_signal, time, verbose=True)
fig.update_yaxes(automargin=True)
fig.update_layout(title=f'powspec @ B = {B:.2f} T, O = {O_string}, T={D_parallel*1000:.0f} MHz',
                  margin=dict(t=75),)
fig.show(renderer='browser')
# fig.write_html(f'Figures/TF_simulations/TF_powder_spectrum_B{B}_O{O_string}_filter{transition_filter}.html')
#%% Plot angular dependence of spectra at fixed B
B = 0.1  # Tesla
thetas = np.linspace(0, 90, 12)
colors = px.colors.qualitative.G10[:len(thetas)]

# Create a common figure
fig = go.Figure()

for theta, color in zip(thetas_deg, colors):
    subfig = stick_spectrum(results=results, theta=theta, B=B)
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
    showlegend=True,
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
# fig.show(renderer='browser')
# fig.write_html('../../Figures/ALC_simulations/TF_B_dependence.html')

#%%
def plot_powder_histogram(results, transition_filter=None, B=0.0, bins=300, normalize=True):
    """
    Plot a rudimentary powder spectrum as a histogram of frequency components.

    Parameters
    ----------
    results : xr.Dataset
        Simulation dataset with variables 'frequencies', 'amplitudes', and 'theta'.
    transition_filter : list[str], optional
        List of transition types to include, if None, all transition types are included.
    B : float
        Magnetic field value to select (in Tesla).
    bins : int
        Number of histogram bins.
    normalize : bool
        Normalize the histogram to its maximum.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        Interactive histogram figure.
    """

    # Select data at given magnetic field
    sel = results.sel(B=B)

    freqs = sel['frequencies'].values
    amps = sel['amplitudes'].values
    ttypes = sel['transition_types'].values
    thetas = np.radians(sel['theta'].values)

    mask = ~np.isnan(freqs)
    freqs = freqs[mask]
    amps = amps[mask]
    ttypes = ttypes[mask]

    if transition_filter is not None:
        mask = np.isin(ttypes, transition_filter)
        freqs, amps, ttypes = freqs[mask], amps[mask], ttypes[mask]

    # sin(theta) weights normalized to sum to 1
    sin_weights = np.sin(thetas)
    sin_weights /= sin_weights.sum()

    # Repeat sin(theta) weights for all transitions per angle
    n_theta, n_trans = sel['frequencies'].shape
    theta_weights = np.repeat(sin_weights, n_trans)

    # Mask same as for freqs
    theta_weights = theta_weights[~np.isnan(sel['frequencies'].values.flatten())]

    # If the filter removed some entries, adjust accordingly
    if len(theta_weights) != len(amps):
        theta_weights = theta_weights[:len(amps)]

    weights = amps * theta_weights

    hist, bin_edges = np.histogram(freqs, bins=bins, weights=weights)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    if normalize and np.max(hist) > 0:
        hist /= np.max(hist)

    # Plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=centers, y=hist, mode='lines', line=dict(width=2)))
    fig.update_layout(
        title=f"Powder histogram spectrum @ B = {B:.2f} T",
        xaxis_title="Frequency / GHz",
        yaxis_title="Relative Intensity",
        template="DemonLab",
        margin=dict(t=75),
    )

    return fig

fig = plot_powder_histogram(results, transition_filter=['ZQ', 'DQ'], B=0.1)
fig.show()

#%% Make pandas df for data exploration
B_selected = 10  # example magnetic field (in Tesla)
sel = results.sel(B=B_selected)

# Extract arrays
freqs = sel["frequencies"].values        # shape: (n_theta, n_transition)
amps = sel["amplitudes"].values
types = sel["transition_types"].values   # same shape, dtype object
thetas = sel["theta"].values

# Collect unique transition type labels
unique_types = np.unique(types[~pd.isna(types)])

# Initialize empty DataFrames
freq_df = pd.DataFrame(index=thetas, columns=unique_types)
amp_df = pd.DataFrame(index=thetas, columns=unique_types)
freq_df.index.name = "theta"
amp_df.index.name = "theta"

# Fill each DataFrame
for i, theta in enumerate(thetas):
    for j in range(freqs.shape[1]):
        ttype = types[i, j]
        if ttype is None or (isinstance(ttype, float) and np.isnan(ttype)):
            continue
        freq_df.loc[theta, ttype] = freqs[i, j]
        amp_df.loc[theta, ttype] = amps[i, j]

#%% Make pandas data frame for energy levels for B dependent energy level plot
theta = 45  # degrees

sel = results['energies'].sel(theta=theta)

energy_df = pd.DataFrame(
    sel.values,
    index=sel['B'].values,
)

# energy_df.to_csv('Python/ALC_simulations/Data/energy_levels_A_iso_1GHzdwdd.csv')

# A_iso and D_parallel are in GHz
fig = px.line(energy_df, x=energy_df.index, y=energy_df.columns, title=f'A={A_iso}, D={D_parallel}, O:{O_string}, θ={theta}°')
fig.update_layout(
    margin=dict(t=45),
)
# fig.write_html(f'Figures/TF_simulations/energy_levels_A_iso_{A_iso}GHz_D_{D_parallel}GHz_theta{theta}.html')
fig.show('browser')


