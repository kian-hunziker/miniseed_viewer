import os
import glob
import numpy as np
from obspy import Stream, Trace, UTCDateTime, read
from typing import Union, Tuple

# absolut path for now
DATA_DIR = '/Users/kianhunziker/Documents/UNI/UNIBAS/MA/seisLM/data/tremor/Tremor_daily'


def get_all_miniseed_files(data_dir: str = DATA_DIR) -> list[str]:
    """
    Get a list of all miniSEED files in the specified data directory.
    The function searches through all subdirectories for files with the .mseed extension.
    Args:
        data_dir (str): Path to the data directory. Defaults to DATA_DIR.

    Returns:
        list[str]: List of miniSEED file paths.
    """
    subdirectories = [d for d in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, d))]
    daily_waveform_files = []
    for subdir in subdirectories:
        daily_waveform_files.extend(glob.glob(os.path.join(data_dir, subdir, '*.mseed')))
    return daily_waveform_files

def load_waveform_multiple_days(
        station: str, 
        waveform_files: list[str], 
        date_span: Union[Tuple[UTCDateTime, UTCDateTime], None] = None,
        fs: int = 100,
        ) -> Stream:
    """
    Load waveform data for a given station from a list of waveform files.
    Merge data from multiple days.
    
    Args:
        station (str): Station code to load data for.
        waveform_files (list[str]): List of waveform file names.
        date_span (tuple[UTCDateTime, UTCDateTime], optional): Start and end times to trim the data. Defaults to None.
    
    Returns:
        Stream: ObsPy Stream object containing the loaded three component waveform data.
    """
    mseed_names = {p.split('/')[-1]: p for p in waveform_files}

    st = Stream()
    for name, file in mseed_names.items():
        if name.startswith('OV.' + station):  # current format: OV.STATION.--.CHANNEL.YYYYDDD.mseed, e.g. OC.INDI.--.HHE.2020237.mseed
            st += read(file)[0]
    if st[0].stats.sampling_rate != fs:
        st = st.resample(fs)
    # optionally: sort, merge, etc.
    # for now, the data seems to be continuous
    st.merge()
    if date_span is not None:
        st.trim(starttime=date_span[0], endtime=date_span[1])

    return st