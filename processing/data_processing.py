import numpy as np
from scipy.signal import hilbert
from obspy import Stream, Trace, UTCDateTime


def clip_trace_amplitude(
        stream: Stream,
        max_amplitude: float = 1e6
):
    """
    Clip the amplitude of each trace in the stream to the specified maximum amplitude.
    """
    for tr in stream:
        tr.data = np.clip(tr.data, -max_amplitude, max_amplitude)
    return stream

def set_high_amplitude_gaps_to_zero(
        stream: Stream,
        amplitude_threshold: float = 25e6
):
    """
    Set data points in each trace of the stream to zero where the absolute amplitude exceeds the specified threshold.
    """
    for tr in stream:
        tr.data = np.where(np.abs(tr.data) > amplitude_threshold, 0, tr.data)
    return stream

# ==============================================================================
# processing functions
# ==============================================================================
def bandpass_filter(stream: Stream, freqmin: float, freqmax: float, fs: int = 100) -> Stream:
    """Apply bandpass filter to each trace in the stream."""
    filtered_stream = stream.copy()
    # taper to reduce edge effects
    # compute taper percentage based on 1s
    taper_percentage = (1.0 * fs) / (len(stream[0].data))
    filtered_stream.taper(max_percentage=taper_percentage, type='cosine')
    filtered_stream.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
    return filtered_stream

def normalize_stream(data: np.ndarray) -> np.ndarray:
    """Normalize the data to have zero mean and unit variance."""
    mean = np.mean(data)
    std = np.std(data)
    std = std if std != 0 else 1.0  # Prevent division by zero
    normalized_data = (data - mean) / std
    return normalized_data

def select_time_window(stream: Stream, starttime: UTCDateTime, endtime: UTCDateTime) -> Stream:
    trimmed_stream = stream.slice(starttime=starttime, endtime=endtime)
    trimmed_stream = trimmed_stream.copy()
    return trimmed_stream

def select_component(stream: Stream, component: str) -> Stream:
    selected_stream = stream.select(component=component)
    try:
        selected_stream = selected_stream.select(channel='HH' + component[-1])
    except Exception:
        print(f'Could not select channel HH{component[-1]}')
        pass
    return selected_stream

def compute_envelope(data: np.ndarray, cutoff=0.1, fs: int = 100) -> np.ndarray:
    analytic_signal = hilbert(data)
    envelope = np.abs(analytic_signal)
    tr = Trace(data=envelope)
    tr.stats.sampling_rate = fs  # Assuming 100 Hz sampling rate
    tr.filter('lowpass', freq=cutoff, corners=4, zerophase=True)
    tr.resample(1)  # downsample to 1 Hz

    return tr.data