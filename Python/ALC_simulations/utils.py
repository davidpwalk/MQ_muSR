import numpy as np
import xarray as xr
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
from scipy.fft import fft, fftfreq, fftshift


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

    if zero_tol is not None:
        small = np.abs(freqs) <= zero_tol
        freqs[small] = 0.0

    adj = np.abs(freqs[:, None] - freqs[None, :]) <= tol

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
    merged_amps = np.array(merged_amps)[order].tolist()
    merged_types = np.array(merged_types, dtype=object)[order].tolist()

    return merged_freqs, merged_amps, merged_types


def calc_powder_signal(results, transition_types=None):
    if transition_types:
        signals = np.zeros_like(results['signal'])
        for transition_type in transition_types:
            if transition_type not in results['transition_type'].values:
                raise ValueError(f"Transition type '{transition_type}' not found in results.")
            else:
                signals += results['signal'].sel(transition_type=transition_type)
        signals = results['signal'].sel(transition_type=transition_type)
    else:
        signals = results['signal'].sum(dim='transition_type')

    weights = np.sin(np.radians(results.theta))
    powder_signals = (signals * weights).sum(dim='theta') / weights.sum()
    return powder_signals


def apodize(signal, time):
    n = len(signal)
    dt = time[1] - time[0]
    window = np.hanning(2*n)[n:]
    apodized_signal = signal * window

    n_fft = 1 << (n - 1).bit_length()  # next power of 2
    padded = np.zeros(n_fft)
    padded[:n] = apodized_signal
    time_padded = np.arange(n_fft) * dt
    return padded, time_padded


def ft(signal, time):
    signal, time = apodize(signal, time)
    fig = px.line(x=time, y=signal, labels={'x': 'Time / ns', 'y': 'Signal'})
    fig.show()
    dt = time[1] - time[0]
    power_spectrum = (np.abs(fft(signal)))**2
    power_spectrum = fftshift(power_spectrum/max(power_spectrum))
    freq = fftshift(fftfreq(len(signal), d=dt))
    return power_spectrum, freq


def plot_powder_spectrum(signal):
    spectrum, freq = ft(signal)
    fig = px.line(x=freq, y=spectrum, labels={'x': 'Frequency / GHz', 'y': 'Intensity'})
    # fig = px.line(x=np.arange(spectrum.shape[0]), y=spectrum, labels={'x': 'Frequency / GHz', 'y': 'Intensity'})
    return fig, freq


# Plotting functions
def time_signal(results, theta, B, transition_type=None):
    if transition_type:
        mask = results['transition_types'].sel(theta=theta, B=B) == transition_type
        freqs = results['frequencies'].sel(theta=theta, B=B).where(mask, drop=True).values
        amps = results['amplitudes'].sel(theta=theta, B=B).where(mask, drop=True).values
    else:
        freqs = results['frequencies'].sel(theta=theta, B=B).values
        amps = results['amplitudes'].sel(theta=theta, B=B).values

    # remove NaN padding
    mask = np.isfinite(freqs) & np.isfinite(amps)
    freqs = freqs[mask]
    amps = amps[mask]

    time = results.coords["time"].values
    signal = np.zeros_like(time, dtype=float)
    for amp, freq in zip(amps, freqs):
        signal += amp * np.cos(2 * np.pi * freq * time)

    fig = px.line(x=time, y=signal, labels={'x': 'Time / ns', 'y': 'Signal'})
    return fig


def stick_spectrum(results, theta, B, transition_type=None, merge_tol=1e-8):
    """
    Display stick spectrum at given θ and B.
    Optionally merges nearby frequencies (within `merge_tol`) only for visualization.
    """
    sel_freqs = results['frequencies'].sel(theta=theta, B=B)
    sel_amps = results['amplitudes'].sel(theta=theta, B=B)
    sel_types = results['transition_types'].sel(theta=theta, B=B)

    # Handle one or multiple transition types
    if transition_type is not None:
        if isinstance(transition_type, (list, tuple, set)):
            mask = xr.DataArray(
                np.isin(sel_types.values, list(transition_type)),
                dims=sel_types.dims,
                coords=sel_types.coords
            )
        else:
            mask = sel_types == transition_type

        freqs = sel_freqs.where(mask, drop=True).values
        amps = sel_amps.where(mask, drop=True).values
        ttypes = sel_types.where(mask, drop=True).values
    else:
        freqs = sel_freqs.values
        amps = sel_amps.values
        ttypes = sel_types.values

    # Remove NaN padding
    mask = np.isfinite(freqs) & np.isfinite(amps)
    freqs, amps, ttypes = freqs[mask], amps[mask], ttypes[mask]

    # ---- Merge nearby lines for cleaner display ----
    if len(freqs) > 0:
        order = np.argsort(freqs)
        freqs, amps, ttypes = freqs[order], amps[order], ttypes[order]

        merged_freqs, merged_amps, merged_types = [], [], []
        current_group = [0]

        for i in range(1, len(freqs)):
            if abs(freqs[i] - freqs[current_group[-1]]) <= merge_tol:
                current_group.append(i)
            else:
                # merge current group
                idx = np.array(current_group)
                total_amp = amps[idx].sum()
                freq_rep = (freqs[idx] * amps[idx]).sum() / total_amp
                merged_freqs.append(freq_rep)
                merged_amps.append(total_amp)
                merged_types.append(",".join(sorted(set(ttypes[idx]))))
                current_group = [i]

        # merge last group
        idx = np.array(current_group)
        total_amp = amps[idx].sum()
        freq_rep = (freqs[idx] * amps[idx]).sum() / total_amp
        merged_freqs.append(freq_rep)
        merged_amps.append(total_amp)
        merged_types.append(",".join(sorted(set(ttypes[idx]))))

        freqs, amps, ttypes = np.array(merged_freqs), np.array(merged_amps), np.array(merged_types)

    # ---- Plot ----
    fig = go.Figure()
    color_cycle = px.colors.qualitative.G10
    for idx, (freq, amp, ttype) in enumerate(zip(freqs, amps, ttypes)):
        color = color_cycle[idx % len(color_cycle)]
        fig.add_trace(
            go.Scatter(
                x=[freq, freq],
                y=[0, amp],
                mode="lines",
                name=f"{freq:.4f} GHz; {amp:.4f}; {ttype}",
                line=dict(color=color),
            )
        )

    fig.add_hline(y=0, line=dict(color="black"), opacity=1)
    fig.update_layout(
        xaxis_title="Frequency / GHz",
        yaxis_title="Amplitude",
        template="DemonLab"
    )
    return fig