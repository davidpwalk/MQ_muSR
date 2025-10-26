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

from Python.jacobi_math import jacobi

# Set custom plotting style
from Python.plot_settings import set_demonlab_style
set_demonlab_style()
pio.templates.default = "DemonLab"

# Import utility functions
from Python.ALC_simulations.utils import *

#%% Parameters and calculations of spin operators and tensors
# Filename and description for saving results (set filename to None to skip saving)
# filename = 'Python/ALC_simulations/Data/TF_NMR_isotropic_cc.nc'
filename = None
desc = 'TF NMR simulation for S=1/2, I=1/2 system with isotropic hyperfine coupling of A_iso=2 MHz, O=Ix'

# Magnetic field range (in Tesla)
magnetic_fields = np.linspace(0, 0.5, 100)
magnetic_fields = [0.1]

# Zeeman (gamma / GHz/T)
ge = 28.02495
gmu = -0.1355

# Coupling constants (in GHz) and rotation angle (in degrees)
# A_iso = 0.5148
A_iso = 0.0
# A_iso = 0.002
D_parallel = 0.002
# D_parallel = 0
D_perp = -D_parallel/2

# Rotation angles (in degrees)
thetas = np.radians(np.linspace(0, 90, 200, dtype=np.float64))
# thetas = np.radians([0, 45, 90, 180])
thetas = [0]

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
O = Ix  # observable
O_string = 'Ix'
threshold = 1e-4  # amplitude threshold for transitions
# TODO: change filter to so all signals are saved and filter only at plotting stage


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

transition_types = ["SQMu", "SQE", "ZQ", "DQ"]

signals_arr = np.zeros((len(transition_types), n_theta, n_B, n_time))

# TODO: fix error where merging transitions leads to error in signal calculation
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
        energies, psi = jacobi(np.real(H_tot.full()))

        energies_arr[ii, jj, :] = energies

        O_dense = O.full()
        O_eigen = psi.conj().T @ O_dense @ psi

        freqs = []
        amps = []
        ttypes = []
        for kk in range(len(energies)):
            for ll in range(kk + 1, len(energies)):
                amp = abs(O_eigen[kk, ll]) ** 2
                ttype = classify_transition(ll + 1, kk + 1)
                # print(f'ttype: {ttype}, kk: {kk+1}, ll: {ll+1}, amp: {amp}, O: \n{O_eigen}') if ttype == 'DQ' or ttype == 'ZQ' else None
                print(f'theta: {np.degrees(thetas[ii])}, energies: {energies}, \n T_lab: {T_lab}')
                if amp > threshold:
                    freqs.append(abs(energies[ll] - energies[kk]))
                    amps.append(amp)
                    ttype = classify_transition(ll + 1, kk + 1)
                    ttypes.append(ttype)

        # store frequencies/amplitudes (pad with NaN if fewer than max_transitions)
        # freqs, amps, ttypes = merge_transitions(freqs, amps, ttypes)
        print(f'freqs: {freqs}, amps: {amps}, types: {ttypes}')
        n_trans = len(freqs)
        frequencies_arr[ii, jj, :n_trans] = freqs
        amplitudes_arr[ii, jj, :n_trans] = amps
        ttypes_arr[ii, jj, :n_trans] = ttypes

        # build signal
        signal_SQMu = np.zeros_like(time, dtype=float)
        signal_SQE = np.zeros_like(time, dtype=float)
        signal_ZQ = np.zeros_like(time, dtype=float)
        signal_DQ = np.zeros_like(time, dtype=float)

        for amp, freq, ttype in zip(amps, freqs, ttypes):
            if ttype == "SQMu":
                signal_SQMu += amp * np.cos(2 * np.pi * abs(freq) * time)
            elif ttype == "SQE":
                signal_SQE += amp * np.cos(2 * np.pi * abs(freq) * time)
            elif ttype == "ZQ":
                signal_ZQ += amp * np.cos(2 * np.pi * abs(freq) * time)
            elif ttype == "DQ":
                signal_DQ += amp * np.cos(2 * np.pi * abs(freq) * time)
            else:
                raise ValueError(f"Unknown transition type: {ttype}")

        type_to_index = {t: i for i, t in enumerate(transition_types)}

        for amp, freq, ttype in zip(amps, freqs, ttypes):
            idx = type_to_index[ttype]
            signals_arr[idx, ii, jj, :] += amp * np.cos(2 * np.pi * abs(freq) * time)

# Wrap into Dataset
results = xr.Dataset(
    {
        "signal": (["transition_type", "theta", "B", "time"], signals_arr),
        "energies": (["theta", "B", "level"], energies_arr),
        "frequencies": (["theta", "B", "transition"], frequencies_arr),
        "amplitudes": (["theta", "B", "transition"], amplitudes_arr),
        "transition_types": (["theta", "B", "transition"], ttypes_arr),
    },
    coords={
        "transition_type": transition_types,
        "theta": np.degrees(thetas),
        "B": magnetic_fields,
        "time": time,
        "level": np.arange(4),
        "transition": np.arange(max_transitions),
    }
)

# TODO: find a smoother solution for saving the data as this can lead to HUGE files if multiple angles and fields are simulated
# TODO: one idea would be to not save time signals and reconstruct them when loading data
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
#%% test
freq_test = 0.002  # GHz
time_test = np.linspace(0, 8000, 16000)  # ns
signal_test = np.cos(2 * np.pi * freq_test * time_test)

px.line(x=time_test, y=signal_test).show()


#%% get signals

def generate_signals(results, time, transition_filter=None):
    """
    Generate time-domain signals from xarray results and return as a pandas DataFrame.

    Parameters
    ----------
    results : xr.Dataset
        Dataset containing 'frequencies', 'amplitudes', and 'transition_types'.
    time : np.ndarray
        Time vector used to reconstruct the signal.
    transition_filter : list[str], optional
        If provided, only include these transition types.

    Returns
    -------
    pd.DataFrame
        Columns correspond to transition types, index corresponds to time.
    """
    # Extract arrays
    freqs = results["frequencies"].values
    amps = results["amplitudes"].values
    ttypes = results["transition_types"].values
    transition_types = results["transition_type"].values

    # Initialize signal accumulator
    signals = {tt: np.zeros_like(time, dtype=float) for tt in transition_types}

    # Iterate over all transitions
    for theta_idx in range(freqs.shape[0]):          # theta dimension
        for B_idx in range(freqs.shape[1]):          # magnetic field dimension
            for tr_idx in range(freqs.shape[2]):     # transition dimension
                freq = freqs[theta_idx, B_idx, tr_idx]
                amp = amps[theta_idx, B_idx, tr_idx]
                ttype = ttypes[theta_idx, B_idx, tr_idx]

                # Skip padded NaN entries
                if np.isnan(freq) or np.isnan(amp) or not isinstance(ttype, str):
                    continue

                if transition_filter and ttype not in transition_filter:
                    continue
                if ttype not in signals:
                    continue

                signals[ttype] += amp * np.cos(2 * np.pi * abs(freq) * time)

    # Convert to DataFrame
    df = pd.DataFrame(signals, index=time)

    signal_total = np.zeros_like(time, dtype=float)
    for ttype in transition_types:
        signal_total += df[ttype]

    df['Total'] = signal_total
    df.index.name = "time"
    return df

signals_df = generate_signals(results, time, transition_filter=['SQE', 'SQMu'])

px.line(x=signals_df.index, y=signals_df['Total']).show()

test = results['signal'].sel(transition_type='SQMu', theta=0, B=0.1)
test += results['signal'].sel(transition_type='SQE', theta=0, B=0.1)

px.line(x=results['time'], y=test).show()

#%% Plotting single spectra
B = 10  # Tesla
theta = 0  # degrees
fig = time_signal(results, theta, B)
fig.update_layout(xaxis_range=[0, 30])
# fig.show()

fig = stick_spectrum(results, theta, B, transition_type=None)
fig.update_layout(
    title=f'TF muSR @ θ = {theta}°, B = {B} T, O = {O_string}',
    margin=dict(t=75),
)
fig.show()
# fig.write_html('../../Figures/ALC_simulations/TF_B10_OIx_D2MHz.html')
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

powder_signals = calc_powder_signal(results, transition_types=None)

fig = px.line(x=powder_signals.sel(B=1, method='nearest').time, y=powder_signals.sel(B=1, method='nearest'), labels={'x': 'Time / ns',
                                                                                                                     'y': 'Signal'})
fig.show()

#%% Plot powder spectra
B = 1  # Tesla
fig, freq = plot_powder_spectrum(powder_signals.sel(B=B, method='nearest').values)
fig.update_yaxes(automargin=True)
fig.update_layout(title=f'Powder muSR spectrum at B = {B} T, O = Ix')
fig.show()
# fig.write_html(f'../../Figures/ALC_simulations/TF_powder_spectrum_thin_pake_pattern_{B}T_Ix_Sx.html')
