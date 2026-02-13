import os
import glob
import numpy as np
from obspy import Stream, Trace, UTCDateTime, read
from typing import Union, Tuple
import pandas as pd
import torch
import json
import re


def get_all_miniseed_files(data_dir: str) -> list[str]:
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
        # TODO: check for valid file extensions
        daily_waveform_files.extend(glob.glob(os.path.join(data_dir, subdir, '*')))
        # find all files (not directories)
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
        file_station = name.split('.')[1]
        if file_station == station:
            st += read(file)[0]
    if len(st) == 0:
        print(f'No data found for station {station} in provided waveform files.')
        # make empty stream with three components and correct sampling rate and start time
        npts = int((date_span[1] - date_span[0]) * fs) if date_span is not None else fs * 60 * 10  # default to 10 minutes of data
        for component in ['HHE', 'HHN', 'HHZ']:
            tr = Trace(data=np.zeros(npts))
            tr.stats.station = station
            tr.stats.channel = component
            tr.stats.sampling_rate = fs
            if date_span is not None:
                tr.stats.starttime = date_span[0]
            st += tr
        
    if st[0].stats.sampling_rate != fs:
        st = st.resample(fs)
    # optionally: sort, merge, etc.
    # for now, the data seems to be continuous
    st.merge(method=1, fill_value=0)
    if date_span is not None:
        st.trim(starttime=date_span[0], endtime=date_span[1], pad=True, fill_value=0)
        st.split().merge(method=1, fill_value=0)

    return st


def read_event_catalog(
        catalog_file: str, 
        start_time_alias: str = 'start_time', 
        end_time_alias: str = 'end_time',
        info_columns: list[str] = ['notes']) -> dict:
    """
    Read an event catalog from a file. Supports JSON, CSV, and PyTorch formats.

    Args:
        catalog_file (str): Path to the catalog file.
        start_time_alias (str, optional): Column name for event start time. Defaults to 'start_time'.
        end_time_alias (str, optional): Column name for event end time. Defaults to 'end_time'.
        info_columns (list[str], optional): List of column names to include as additional info. Defaults to ['notes'].
    Returns:
        dict: Dictionary of events with start time, end time, and additional info.
    """

    file_extension = catalog_file.split('.')[-1].lower()
    if file_extension == 'json':
        with open(catalog_file, 'r') as f:
            catalog = json.load(f)
    elif file_extension == 'csv':
        df = pd.read_csv(catalog_file)
        # add index as first column if not present
        if df.columns[0] != 'event_id':
            df.insert(0, 'event_id', range(len(df)))
        catalog = df.to_dict(orient='index')
    elif file_extension in ['pt', 'pth']:
        catalog = torch.load(catalog_file)
    else:
        raise ValueError(f"Unsupported catalog file format: {file_extension}")
    
    print(f"Loaded {len(catalog)} events from catalog at {catalog_file}")
    """
    print(catalog.keys())
    for k, v in catalog.items():
        print(f"Event {k}: {v}")
    """
    output_catalog = {}
    for i, event in enumerate(catalog.values()):
        start_time = event[start_time_alias]
        end_time = event[end_time_alias]
        if isinstance(start_time, str):
            start_time = UTCDateTime(start_time)
        if isinstance(end_time, str):
            end_time = UTCDateTime(end_time)
        notes = ''
        for col in info_columns:
            notes += str(event.get(col, '')) + ', '
        notes = notes.rstrip(', ')
        output_catalog[i] = {
            'start_time': start_time,
            'end_time': end_time,
            'notes': notes
        }
        # add remaining columns as additional info
        for col in event.keys():
            if col not in [start_time_alias, end_time_alias] + info_columns:
                output_catalog[i][col] = event[col]
    return output_catalog

def write_event_catalog(
        event_dict: dict,
        output_file: str
):
    """
    Write event information to a file. Supports JSON, CSV, and PyTorch formats.

    Args:
        event_dict (dict): Dictionary of events to write.
        output_file (str): Path to the output JSON file.
    """
    serializable_dict = {}
    for event_id, event in event_dict.items():
        serializable_dict[event_id] = {
            'start_time': str(event['start_time']),
            'end_time': str(event['end_time']),
            'notes': event['notes'].replace('\n', ' ').replace(',', ';')
        }
        # add remaining info
        for key, value in event.items():
            if key not in ['start_time', 'end_time', 'notes']:
                serializable_dict[event_id][key] = value
    file_extension = output_file.split('.')[-1].lower()
    if file_extension == 'json':
        with open(output_file, 'w') as f:
            json.dump(serializable_dict, f, indent=4)
    elif file_extension in ['pt', 'pth']:
        torch.save(serializable_dict, output_file)
    elif file_extension == 'csv':
        df = pd.DataFrame.from_dict(serializable_dict, orient='index')
        df.to_csv(output_file, index_label='event_id')
    else:
        raise ValueError(f"Unsupported output file format: {file_extension}")
    
    print(f"Wrote {len(event_dict)} events to catalog at {output_file}")
    

def date_to_jul_string(date: UTCDateTime) -> str:
    """
    Convert a UTCDateTime object to a Julian date string in the format YYYYDDD.
    
    Args:
        date (UTCDateTime): The date to convert.
    Returns:
        str: The Julian date string.
    """
    year = date.year
    day_of_year = date.julday
    jul_string = f"{year}{day_of_year:03d}"
    return jul_string

def get_jul_strings_for_date_range(start_date: UTCDateTime, end_date: UTCDateTime) -> list[str]:
    """
    Get a list of Julian date strings for each day in the specified date range.
    
    Args:
        start_date (UTCDateTime): The start date.
        end_date (UTCDateTime): The end date.
    Returns:
        list[str]: List of Julian date strings in the format YYYYDDD.
    """
    jul_strings = []
    current_date = start_date
    while current_date <= end_date:
        jul_strings.append(date_to_jul_string(current_date))
        current_date += 86400  # increment by one day in seconds
    return jul_strings

def get_filenames_containing_jul_strings(
        jul_strings: list[str], 
        all_filenames: list[str]
    ) -> list[str]:
    """
    Get a list of filenames that contain any of the specified Julian date strings.
    
    Args:
        jul_strings (list[str]): List of Julian date strings to search for.
        all_filenames (list[str]): List of all filenames to search within.
    
    Returns:
        list[str]: List of filenames that contain any of the Julian date strings.
    """
    matching_filenames = []
    for filename in all_filenames:
        jst = filename.split('.')[-2]
        if jst in jul_strings:
            matching_filenames.append(filename)
    return matching_filenames