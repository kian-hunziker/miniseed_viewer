import sys
import numpy as np
from PySide6 import QtWidgets
from PySide6 import QtCore
import pyqtgraph as pg
import os
import glob
import numpy as np
from obspy import UTCDateTime, Stream, Trace
from pyqtgraph.exporters import ImageExporter, SVGExporter
from tqdm import tqdm
import json
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QFileDialog,
    QMessageBox,
)
from PySide6.QtGui import QAction

from processing.data_loading import (
    get_all_miniseed_files, 
    load_waveform_multiple_days,
    read_event_catalog,
    write_event_catalog,
    date_to_jul_string,
    get_jul_strings_for_date_range,
    get_filenames_containing_jul_strings,
    )

from plotting.mutlistation_plots import MultiStationSingleComponentView
from plotting.zne_plots import ThreeComponentSingleStationView, ThreeComponentMotionPlotView
from plotting.spectogram_plots import MultistationWithSpectogram

from plotting.spectrogram_control import SpectrogramControlWindow

from processing.data_processing import (
    clip_trace_amplitude,
    set_high_amplitude_gaps_to_zero,
)

from plotting.utility_classes import CatalogDialog, EventNoteDialog, CatalogEditor

# ==============================================================================
# MAIN APPLICATION
# ==============================================================================

class TremorViewer(QtWidgets.QMainWindow):
    def __init__(self,
                 data_dir,
                 stations: list[str] = ['JUDI', 'LAFE', 'JACO', 'INDI', ],
                 fs: int = 100,
                 max_npts: int = 720_000,
                 clip_amplitude: float = 20e6,
                 zero_amplitude: float = 25e6,
                 init_from_state: str = 'default', # one of default, last
                 ):
        super().__init__()
        self.data_dir = data_dir
        self.fs = fs
        self.stations = stations
        self.max_npts = max_npts
        self.default_y_range = 4
        self.zero_amplitude = zero_amplitude
        self.clip_amplitude = clip_amplitude

        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        
        # load waveform data for all stations
        self.waveform_files = get_all_miniseed_files(data_dir=self.data_dir)
        self.station_streams = {}
        '''
        for station in tqdm(self.stations, desc="Loading stations"):
            st = load_waveform_multiple_days(
                station=station,
                waveform_files=self.waveform_files,
                fs=self.fs,
            )
            self.station_streams[station] = st

        # clip amplitude
        # also set high amplitude gaps to zero
        for station in self.stations:
            st = self.station_streams[station]
            st = set_high_amplitude_gaps_to_zero(st, amplitude_threshold=self.zero_amplitude)
            st = clip_trace_amplitude(st, max_amplitude=self.clip_amplitude)
            self.station_streams[station] = st'''
        
        # make folder structure
        '''
        -config
        -data
            -events
            -temp_files
        -outputs
        '''
        os.makedirs('config', exist_ok=True)
        os.makedirs('data/events', exist_ok=True)
        os.makedirs('data/temp_files', exist_ok=True)
        os.makedirs('outputs', exist_ok=True)
        
        # setup window
        self.setWindowTitle("Tremor Viewer v3.0")
        self.resize(1200, 600)

        # initial state
        if init_from_state in ['default', 'last']:
            try:
                self._load_state(name=init_from_state)
                #print(self.state)
            except:
                print(f"Failed to load state from config/{init_from_state}.json, using default state.")
                self._setup_initial_state()
        else:
            self._setup_initial_state()

        # load catalogue
        cat_path = self.state.get('last_catalog_load_folder', None)
        cat_file = self.state.get('last_loaded_catalog_file', None)
        if self.state.get('all_catalog_changes_saved', True) == False:
            # if there are unsaved changes
            reply = QMessageBox.question(
                self,
                "Unsaved Catalog Changes",
                "There are unsaved changes to the event catalog. Do you want to restore the last catalog you were working on?",
                QMessageBox.Yes | QMessageBox.No,
            )
            # if yes, restore from temp file
            temp_path = os.path.join('data/temp_files', 'temp_catalog.csv')
            if reply == QMessageBox.Yes and os.path.exists(temp_path):
                    self.event_catalog = read_event_catalog(temp_path)
                    print(f"Restored event catalog from last session: {temp_path}")
            else:
                # try to load from last loaded catalog if temp file doesn't exist
                if cat_file and os.path.exists(os.path.join(cat_path, cat_file)):
                    self.event_catalog = read_event_catalog(os.path.join(cat_path, cat_file))
                    print(f"Loaded previous event catalog from: {os.path.join(cat_path, cat_file)}")
                else:
                    self.event_catalog = {}

        elif cat_file and os.path.exists(os.path.join(cat_path, cat_file)):
            self.event_catalog = read_event_catalog(os.path.join(cat_path, cat_file))
            print(f"Loaded previous event catalog from: {os.path.join(cat_path, cat_file)}")
        else:
            try:
                self.event_catalog_path = 'data/events/tremor_multi_tremor_bin_class_seisLM_base_2025-12-05-12h-41m-26s.csv'
                self.event_catalog = read_event_catalog(self.event_catalog_path)
            except Exception as e:
                print(f"Failed to load event catalog from {self.event_catalog_path}: {e}")
                self.event_catalog = {}
        
        self.load_waveform_data()
            
        self.wf_color = 'k'#(100, 100, 100)  # waveform color
        self.wf_linewidth = 0.8  # waveform line width
        self.env_color = 'r'  # envelope color
        self.env_linewidth = 2  # envelope line width
        self.background_color = 'w' #(180, 180, 180)  # dark background
        
        pg.setConfigOptions(background=self.background_color, foreground=self.wf_color)

        plot_args = {
            'wf_color': self.wf_color,
            'wf_linewidth': self.wf_linewidth,
            'env_color': self.env_color,
            'env_linewidth': self.env_linewidth,
        }

        self.layout = QtWidgets.QVBoxLayout(central_widget)

        # stacked views
        self.stack = QtWidgets.QStackedWidget()

        self.multi_station_view = MultiStationSingleComponentView(
                stations=self.stations,
                state=self.state, 
                stream_dict=self.station_streams,
                max_npts=self.max_npts,
                fs=self.fs,
                **plot_args
            )
        self.multi_station_view_spectrograms = MultistationWithSpectogram(
            stations=self.stations,
            state=self.state, 
            stream_dict=self.station_streams,
            max_npts=self.max_npts,
            fs=self.fs,
            **plot_args
        )
        self.three_comp_view = ThreeComponentSingleStationView(
            stations=self.stations,
            state=self.state,
            stream_dict=self.station_streams,
            max_npts=self.max_npts,
            fs=self.fs,
            **plot_args
            )
        self.three_comp_with_motion_view = ThreeComponentMotionPlotView(
            stations=self.stations,
            state=self.state, 
            stream_dict=self.station_streams,
            max_npts=self.max_npts,
            fs=self.fs,
            **plot_args
            )
        self.stack.addWidget(self.multi_station_view)
        self.stack.addWidget(self.three_comp_view)
        self.stack.addWidget(self.three_comp_with_motion_view)
        self.stack.addWidget(self.multi_station_view_spectrograms)
        self.all_plots = self.multi_station_view.plots + self.three_comp_view.plots + self.three_comp_with_motion_view.plots + self.multi_station_view_spectrograms.plots + self.multi_station_view_spectrograms.spectrogram_plots

        # default view
        if self.state['current_view'] == 'three_comp' and not self.state['show_motion_plot']:
            self.stack.setCurrentWidget(self.three_comp_view)
        elif self.state['current_view'] == 'three_comp_with_motion' and self.state['show_motion_plot']:
            self.stack.setCurrentWidget(self.three_comp_with_motion_view)
        elif self.state['current_view'] == 'multi_station':
            self.stack.setCurrentWidget(self.multi_station_view)
        elif self.state['current_view'] == 'multi_station_spectrograms':
            self.stack.setCurrentWidget(self.multi_station_view_spectrograms)
        else:
            raise ValueError(f"Invalid current_view: {self.state['current_view']}")
        self.layout.addWidget(self.stack)
        
        # build controls
        self._build_controls(self.layout)
        self._create_menu()
        
        # initial refresh
        QtCore.QTimer.singleShot(0, self.stack.currentWidget().refresh)
        
        # setup dragging signals
        self._setup_dragging_signals()
        # setup mouse signals
        self._setup_mouse_signals()

        self.spec_controls = SpectrogramControlWindow(self.state)
        self.spec_controls.paramsChanged.connect(self._update_plots)

        if self.state['spectrogram']['show_control']:
            self.spec_controls.show()

    def _create_menu(self):
        menu_bar = self.menuBar()

        # ---- Data Menu ----
        data_menu = menu_bar.addMenu("Data")

        set_path_action = QAction("Set Data Path…", self)
        set_path_action.triggered.connect(self.select_data_path)
        set_path_action.setStatusTip("Select the folder containing your miniseed data files")
        data_menu.addAction(set_path_action)

        reload_files_action = QAction("Reload Data Files", self)
        reload_files_action.triggered.connect(self.reload_data_files)
        reload_files_action.setStatusTip("Reload the list of miniseed files from the current data folder")
        data_menu.addAction(reload_files_action)

    def select_data_path(self):
        dialog = QFileDialog(self)
        dialog.setWindowTitle("Select Data Folder")
        dialog.setFileMode(QFileDialog.Directory)
        dialog.setOption(QFileDialog.ShowDirsOnly, True)

        # Start in previous path if available
        if self.data_dir and os.path.exists(self.data_dir):
            dialog.setDirectory(self.data_dir)

        if dialog.exec():  # User clicked "Open"
            selected_folder = dialog.selectedFiles()[0]

            # Confirm selection
            reply = QMessageBox.question(
                self,
                "Confirm Data Path",
                f"Use this folder as data path?\n\n{selected_folder}",
                QMessageBox.Yes | QMessageBox.Cancel,
            )

            if reply == QMessageBox.Yes:
                self.data_dir = selected_folder
                self.waveform_files = get_all_miniseed_files(data_dir=self.data_dir)
                print(f"New data path set to: {self.data_dir}")
        else:
            print("Selection canceled")
    
    def reload_data_files(self):
        if self.data_dir and os.path.exists(self.data_dir):
            self.waveform_files = get_all_miniseed_files(data_dir=self.data_dir)
            print(f"Reloaded data files from: {self.data_dir}")
        else:
            print("No valid data directory set. Please set the data path first.")
        
    def _setup_initial_state(self):
        self.state = {
                'component': 'N',
                'start_time': UTCDateTime(2020, 8, 24, 5, 0),
                'end_time': UTCDateTime(2020, 8, 24, 7, 0),
                'current_date': UTCDateTime(2020, 8, 24),
                'freqmin': None,
                'freqmax': None,
                'y_limits': {f'{s}_{c}': (-self.default_y_range, self.default_y_range) for s in self.stations for c in ['Z', 'N', 'E']},
                'show_envelope': False,
                'show_waveform': True,
                'show_motion_plot': True,
                'env_cutoff': 0.1,
                'zne_station': self.stations[0],
                'selection_start': None,
                'selection_end': None,
                'dragging': False,
                'current_view': 'multi_station',
                'downsampling_factor': 1,
                'last_plot_save_folder': None,
                'last_catalog_save_folder': None,
                'last_catalog_load_folder': None,
                'spectrogram': {
                    'nperseg': 256,
                    'overlap': 0.5,
                    'cmap': 'viridis',
                    'mode': 'psd',
                    'log_scale': True,
                    'show_control': True,
                    'independent_filtering': False,
                    'f_bandpass_min': 0.1,
                    'f_bandpass_max': 20.0,
                }
            }
        # override saved default state, to update defaults
        self._save_state(name='default')
        
    def load_waveform_data(self):
        """
        check if start_time and end_time are in loaded data
        load new / more data if necessary
        """
        current_start_date = UTCDateTime(self.state['start_time'].year, self.state['start_time'].month, self.state['start_time'].day)
        current_end_date = UTCDateTime(self.state['end_time'].year, self.state['end_time'].month, self.state['end_time'].day)

        has_loaded_data = len(self.station_streams) > 0

        if has_loaded_data:
            data_start_date = min([st[0].stats.starttime.date for st in self.station_streams.values()])
            data_end_date = max([st[0].stats.endtime.date for st in self.station_streams.values()])

        if not has_loaded_data or current_start_date.date < data_start_date or current_end_date.date > data_end_date:
            #jul_strings = get_jul_strings_for_date_range(current_start_date, current_end_date)
            #needed_files = get_filenames_containing_jul_strings(jul_strings=jul_strings, all_filenames=self.waveform_files)
            
            needed_files = []
            day = current_start_date
            while day <= current_end_date + 24*3600:  # add one day to include end date
                # format day into sting 'YYYY-MM-DD'
                day_str = UTCDateTime(day).strftime('%Y-%m-%d')
                # find files whos parent folder match this string
                day_files = [f for f in self.waveform_files if day_str in f]
                needed_files.extend(day_files)
                day += 24*3600
            
            start_time_stamp = UTCDateTime(current_start_date.year, current_start_date.month, current_start_date.day, 0, 0)
            end_time_stamp = UTCDateTime(current_end_date.year, current_end_date.month, current_end_date.day, 23, 59, 59)
            
            print(f'Loading data for: {start_time_stamp} to {end_time_stamp}')

            for station in tqdm(self.stations, desc="Loading additional data for stations"):
                st = load_waveform_multiple_days(
                    station=station,
                    waveform_files=needed_files,
                    date_span=(start_time_stamp, end_time_stamp),
                    fs=self.fs,
                )
                self.station_streams[station] = st
            for station in self.stations:
                st = self.station_streams[station]
                st = set_high_amplitude_gaps_to_zero(st, amplitude_threshold=self.zero_amplitude)
                st = clip_trace_amplitude(st, max_amplitude=self.clip_amplitude)
                self.station_streams[station] = st
            


    def _get_all_plots(self):
        return self.all_plots
    
    def _setup_dragging_signals(self):
        for p in self._get_all_plots():
            vb = p.getViewBox()
            vb.dragging_signal.connect(self.on_dragging_changed)
    
    def _setup_mouse_signals(self):
        # All plots share the same GraphicsScene
        for p in self._get_all_plots():
            vb = p.getViewBox()
            vb.scene().sigMouseClicked.connect(self._on_scene_mouse_clicked)
        scene = self.all_plots[0].scene()
        scene.sigMouseClicked.connect(self._on_scene_mouse_clicked)
    
    def _build_controls(self, layout):
        controls = QtWidgets.QGridLayout()
        layout.addLayout(controls)

        # Time controls
        controls.addWidget(QtWidgets.QLabel("Start (HH:MM or ISO):"), 0, 0)
        self.box_start = QtWidgets.QLineEdit(self.state["start_time"].strftime("%H:%M"))
        controls.addWidget(self.box_start, 0, 1)

        controls.addWidget(QtWidgets.QLabel("End (HH:MM or ISO):"), 0, 2)
        self.box_end = QtWidgets.QLineEdit(self.state["end_time"].strftime("%H:%M"))
        controls.addWidget(self.box_end, 0, 3)

        self.btn_apply_time = QtWidgets.QPushButton("Apply Time")
        controls.addWidget(self.btn_apply_time, 0, 4)

        # Date controls
        controls.addWidget(QtWidgets.QLabel("Date (MM-DD or YYYY-MM-DD):"), 1, 0)
        self.box_date = QtWidgets.QLineEdit(self.state["current_date"].strftime("%m-%d"))
        controls.addWidget(self.box_date, 1, 1)
        self.btn_apply_date = QtWidgets.QPushButton("Apply Date")
        controls.addWidget(self.btn_apply_date, 1, 2)

        # Bandpass controls
        f_min = self.state['freqmin'] if self.state['freqmin'] is not None else 2
        controls.addWidget(QtWidgets.QLabel("Fmin (Hz):"), 2, 0)
        self.box_fmin = QtWidgets.QLineEdit(str(f_min))
        controls.addWidget(self.box_fmin, 2, 1)
        controls.addWidget(QtWidgets.QLabel("Fmax (Hz):"), 2, 2)
        f_max = self.state['freqmax'] if self.state['freqmax'] is not None else 8
        self.box_fmax = QtWidgets.QLineEdit(str(f_max))
        controls.addWidget(self.box_fmax, 2, 3)
        self.btn_apply_bp = QtWidgets.QPushButton("Apply Bandpass")
        controls.addWidget(self.btn_apply_bp, 2, 4)

        # Component radio buttons
        comp_layout = QtWidgets.QHBoxLayout()
        controls.addLayout(comp_layout, 3, 0, 1, 5)
        comp_layout.addWidget(QtWidgets.QLabel("Component:"))
        self.rbZ = QtWidgets.QRadioButton("Z")
        self.rbN = QtWidgets.QRadioButton("N")
        self.rbE = QtWidgets.QRadioButton("E")
        self.rbZNE = QtWidgets.QRadioButton("ZNE View")
        comp_layout.addWidget(self.rbZ)
        comp_layout.addWidget(self.rbN)
        comp_layout.addWidget(self.rbE)
        comp_layout.addWidget(self.rbZNE)
        
        if self.state['current_view'] == 'multi_station' or self.state['current_view'] == 'multi_station_spectrograms':
            self.rbZNE.setChecked(False)
            self.rbZ.setChecked(self.state['component'] == 'Z')
            self.rbN.setChecked(self.state['component'] == 'N')
            self.rbE.setChecked(self.state['component'] == 'E')
        else:
            self.rbZNE.setChecked(True)
        
        # Station selector for ZNE view
        self.station_selector = QtWidgets.QComboBox()
        self.station_selector.addItems(self.stations)  # populate with all station names
        comp_layout.addWidget(self.station_selector)
        self.station_selector.setCurrentText(self.state['zne_station'])

        # event selector
        self.event_list = QtWidgets.QComboBox()
        event_names = [
                self._parse_event_name(e, i)
                for i, e in self.event_catalog.items()
            ]
        self.event_list.addItems(event_names)
        controls.addWidget(self.event_list, 4, 0, 1, 2)

        # envelope and waveform visibility
        self.cb_show_spectrogram = QtWidgets.QCheckBox("Show Spectrogram")
        self.cb_show_spectrogram.setChecked(self.state['current_view'] == 'multi_station_spectrograms')
        controls.addWidget(self.cb_show_spectrogram, 0, 6)
        self.cb_show_waveform = QtWidgets.QCheckBox("Show Waveform")
        self.cb_show_waveform.setChecked(self.state['show_waveform'])
        controls.addWidget(self.cb_show_waveform, 0, 5)
        self.cb_show_envelope = QtWidgets.QCheckBox("Show Envelope")
        self.cb_show_envelope.setChecked(self.state['show_envelope'])
        controls.addWidget(self.cb_show_envelope, 1, 5)
        self.cb_show_motion_plot = QtWidgets.QCheckBox("Show Motion Plot")
        self.cb_show_motion_plot.setChecked(self.state['show_motion_plot'])
        controls.addWidget(self.cb_show_motion_plot, 0, 7)
        self.show_spectrogram_control = QtWidgets.QCheckBox("Spectrogram Settings")
        controls.addWidget(self.show_spectrogram_control, 2, 5)
        self.show_spectrogram_control.setChecked(self.state['spectrogram']['show_control'])
        
        self.box_env_cutoff = QtWidgets.QLineEdit(str(self.state['env_cutoff']))
        controls.addWidget(self.box_env_cutoff, 1, 6)
        self.btn_apply_env_cutoff = QtWidgets.QPushButton("Apply Env Cutoff")
        controls.addWidget(self.btn_apply_env_cutoff, 1, 7)

        # Reset Y button
        self.btn_reset_y = QtWidgets.QPushButton("Reset Y-Limits")
        controls.addWidget(self.btn_reset_y, 4, 5)

        # save button
        self.btn_save = QtWidgets.QPushButton("Save Figure")
        controls.addWidget(self.btn_save, 4, 7)

        # add event button
        self.btn_add_event = QtWidgets.QPushButton("Add Event")
        controls.addWidget(self.btn_add_event, 4, 2)
        self.btn_save_load_catalog = QtWidgets.QPushButton("Save/Load Catalog")
        controls.addWidget(self.btn_save_load_catalog, 4, 3)
        self.btn_edit_catalog = QtWidgets.QPushButton("Edit Catalog")
        controls.addWidget(self.btn_edit_catalog, 4, 4)

        # Connect signals
        self.station_selector.currentTextChanged.connect(self.on_station_changed)
        self.event_list.currentTextChanged.connect(self.on_event_changed)
        self.btn_apply_bp.clicked.connect(self.on_apply_bandpass)
        self.btn_apply_time.clicked.connect(self.on_apply_time)
        self.btn_apply_date.clicked.connect(self.on_apply_date)
        self.rbZ.toggled.connect(self.on_component_changed)
        self.rbN.toggled.connect(self.on_component_changed)
        self.rbE.toggled.connect(self.on_component_changed)
        self.rbZNE.toggled.connect(self.change_view_mode)
        self.btn_reset_y.clicked.connect(self.on_reset_y)
        self.btn_save.clicked.connect(self.save_plot)
        self.cb_show_envelope.stateChanged.connect(self.show_envelope_changed)
        self.cb_show_waveform.stateChanged.connect(self.show_waveform_changed)
        self.cb_show_motion_plot.stateChanged.connect(self.show_motion_changed)
        self.btn_apply_env_cutoff.clicked.connect(self.on_apply_env_cutoff)
        self.btn_add_event.clicked.connect(self.on_event_added)
        self.btn_save_load_catalog.clicked.connect(self.catalog_dialog)
        self.btn_edit_catalog.clicked.connect(self.edit_catalog)
        self.cb_show_spectrogram.stateChanged.connect(self.change_view_mode)
        self.show_spectrogram_control.stateChanged.connect(self.on_show_spectrogram_control)

        # connect return pressed in text boxes to apply functions
        self.box_start.returnPressed.connect(self.on_apply_time)
        self.box_end.returnPressed.connect(self.on_apply_time)
        self.box_date.returnPressed.connect(self.on_apply_date)
        self.box_fmin.returnPressed.connect(self.on_apply_bandpass)
        self.box_fmax.returnPressed.connect(self.on_apply_bandpass)
        self.box_env_cutoff.returnPressed.connect(self.on_apply_env_cutoff)
    
    def _parse_time_input(self, text):
        text = text.strip()
        if not text:
            return None
        # try ISO first
        try:
            time = UTCDateTime(text)
            return time
        except Exception:
            pass
        # try HH:MM or HH:MM:SS using base_date
        try:
            parts = text.split(":")
            if len(parts) >= 2:
                hh = int(parts[0])
                mm = int(parts[1])
                ss = int(parts[2]) if len(parts) >= 3 else 0
                base = self.state["current_date"]
                return base + hh * 3600 + mm * 60 + ss
        except Exception:
            return None
        return None
    
    def _parse_date_string(self, s):
        parts = s.split("-")
        if len(parts) == 2:
            month = int(parts[0])
            day = int(parts[1])
            # default year kept as 2020
            return UTCDateTime(f"2020-{month:02d}-{day:02d}")
        elif len(parts) == 3:
            year = int(parts[0])
            month = int(parts[1])
            day = int(parts[2])
            return UTCDateTime(f"{year}-{month:02d}-{day:02d}")
        else:
            raise ValueError("Date must be MM-DD or YYYY-MM-DD")
        
    def _save_state(self, name='last'):
        """
        Safe state as json file to config/name.json
        """
        state_copy = self.state.copy()
        # convert UTCDateTime to isoformat strings
        for key in state_copy.keys():
            if state_copy[key] is not None:
                if isinstance(state_copy[key], UTCDateTime):
                    state_copy[key] = state_copy[key].isoformat()
        with open(f"config/{name}.json", "w") as f:
            json.dump(state_copy, f, indent=4)

    def _load_state(self, name='default'):
        """
        Load state from config/state.json
        """
        with open(f"config/{name}.json", "r") as f:
            state_loaded = json.load(f)
        utc_keys = ['start_time', 'end_time', 'current_date', 'selection_start', 'selection_end']
        for key in state_loaded.keys():
            if state_loaded[key] is not None:
                if key in utc_keys and isinstance(state_loaded[key], str):
                    state_loaded[key] = UTCDateTime(state_loaded[key])
        self.state = state_loaded

    def _update_plots(self):
        self._save_state(name='last')
        self.load_waveform_data()
        current_view = self.stack.currentWidget()
        current_view.refresh()

    def on_dragging_changed(self, dragging):
        # update dragging state
        self.state['dragging'] = dragging

    def on_show_spectrogram_control(self, state):
        if state > 0:
            self.spec_controls.show()
        else:
            self.spec_controls.hide()

    def change_view_mode(self):
        if self.rbZNE.isChecked():
            print("Switching to 3-component view")
            if self.state['show_motion_plot']:
                self.state['current_view'] = 'three_comp_with_motion'
                self.stack.setCurrentWidget(self.three_comp_with_motion_view)
            else:
                self.state['current_view'] = 'three_comp'
                self.stack.setCurrentWidget(self.three_comp_view)
        else:
            if self.cb_show_spectrogram.isChecked():
                print("Switching to multi-station spectrogram view")
                self.state['current_view'] = 'multi_station_spectrograms'
                self.stack.setCurrentWidget(self.multi_station_view_spectrograms)
            else:
                print("Switching to multi-station single-component view")
                self.state['current_view'] = 'multi_station'
                self.stack.setCurrentWidget(self.multi_station_view)
        self._update_plots()
    
    def on_apply_bandpass(self):
        try:
            fmin = float(self.box_fmin.text()) if self.box_fmin.text().strip() != "" else None
            fmax = float(self.box_fmax.text()) if self.box_fmax.text().strip() != "" else None
            self.state['freqmin'] = fmin
            self.state['freqmax'] = fmax
            print(f"Applying bandpass filter: {fmin} - {fmax} Hz")
            self._update_plots()
        except ValueError:
            pass  # ignore invalid input
    
    def on_reset_y(self):
        # loop through y-limits and reset to default
        current_widget = self.stack.currentWidget()
        for p in current_widget.plots:
            p.setYRange(-self.default_y_range, self.default_y_range)
        self._update_plots()

    def on_apply_time(self):
        t0 = self._parse_time_input(self.box_start.text())
        t1 = self._parse_time_input(self.box_end.text())
        if t0 and t1 and t1 > t0:
            self.state['start_time'] = t0
            self.state['end_time'] = t1
            self._update_plots()

    def on_component_changed(self):
        if self.rbZ.isChecked():
            self.state['component'] = 'Z'
        elif self.rbN.isChecked():
            self.state['component'] = 'N'
        else:
            self.state['component'] = 'E'
        self._update_plots()

    def show_envelope_changed(self, state):
        """Toggle envelope visibility."""
        print(f"Show envelope: {state > 0}")
        self.state['show_envelope'] = state > 0
        self._update_plots()
    
    def show_waveform_changed(self, state):
        """Toggle waveform visibility."""
        print(f"Show waveform: {state > 0}")
        self.state['show_waveform'] = state > 0
        self._update_plots()

    def show_motion_changed(self, state):
        "toggle motion plot visibility in 3-component view"
        self.state['show_motion_plot'] = state > 0
        if self.state['current_view'] in ['three_comp', 'three_comp_with_motion']:
            print(f"show motion plot: {state > 0}")
            if state > 0:
                self.stack.setCurrentWidget(self.three_comp_with_motion_view)
                self.state['current_view'] = 'three_comp_with_motion'
            else:
                self.stack.setCurrentWidget(self.three_comp_view)
                self.state['current_view'] = 'three_comp'
            self._update_plots()
    
    def on_apply_date(self):
        date = self._parse_date_string(self.box_date.text())
        self.state["current_date"] = date
        # update start/end times to new date but same HH:MM:SS
        start_time = self.state["start_time"]
        end_time = self.state["end_time"]
        new_start = UTCDateTime(date.year, date.month, date.day, start_time.hour, start_time.minute, start_time.second)
        new_end = UTCDateTime(date.year, date.month, date.day, end_time.hour, end_time.minute, end_time.second)
        self.state["start_time"] = new_start
        self.state["end_time"] = new_end
        self._update_plots()
    
    def mouseReleaseEvent(self, event):
        # get current x-limits from first plot
        if self.state['dragging']:
            self.state['dragging'] = False
            current_view = self.stack.currentWidget()
            if current_view.shift_pressed:
                return  # ignore if shift is pressed
            
            x_min, x_max = self.stack.currentWidget().plots[0].viewRange()[0]
            start_time = self.state['start_time'] + self.state['downsampling_factor'] * x_min
            end_time = self.state['start_time'] + self.state['downsampling_factor'] * x_max
            self.state['start_time'] = start_time
            self.state['end_time'] = end_time
            self.state['current_date'] = UTCDateTime(start_time.year, start_time.month, start_time.day)
        
        else:
            new_start = self.state['selection_start']
            new_end = self.state['selection_end']
            # check for selection regions in all plots
            for p in self._get_all_plots():
                vb = p.getViewBox()
                if getattr(vb, 'lr', None) is not None:
                    region = vb.lr.getRegion()
                    x0, x1 = region
                    if x0 != self.state['selection_start']:
                        new_start = self.state['start_time'] + self.state['downsampling_factor'] * x0
                    if x1 != self.state['selection_end']:
                        new_end = self.state['start_time'] + self.state['downsampling_factor'] * x1

                vb.clear_selection()
            self.state['selection_start'] = new_start    
            self.state['selection_end'] = new_end
        
        self._update_plots()
        super().mouseReleaseEvent(event)
    

    def save_plot(self):
        default_name = f'tremor_plot_{self.state["current_date"].strftime("%Y%m%d")}_{self.state["component"]}_{self.state["freqmin"]}-{self.state["freqmax"]}Hz'
        if self.state['last_plot_save_folder']:
            start_path = os.path.join(self.state['last_plot_save_folder'], default_name)
        else:
            start_path = os.path.join("outputs", default_name)
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            None, "Save plot", start_path, "PNG (*.png);;SVG (*.svg)"
        )
        # update last save folder
        if path:
            self.state['last_plot_save_folder'] = os.path.dirname(path)
        if not path:
            return

        current_view = self.stack.currentWidget()

        if path.endswith(".png"):
            ImageExporter(current_view.graphics.scene()).export(path)
            print(f"Saving plot to {path}.png")
        elif path.endswith(".svg"):
            SVGExporter(current_view.graphics.scene()).export(path)
            print(f"Saving plot to {path}.svg")
        # update state with last save folder
        self._save_state(name='last')

    def on_apply_bandpass(self):
        try:
            fmin = float(self.box_fmin.text()) if self.box_fmin.text().strip() != "" else None
            fmax = float(self.box_fmax.text()) if self.box_fmax.text().strip() != "" else None
            self.state['freqmin'] = fmin
            self.state['freqmax'] = fmax
            print(f"Applying bandpass filter: {fmin} - {fmax} Hz")
            self._update_plots()
        except ValueError:
            pass  # ignore invalid input
    
    def on_apply_env_cutoff(self):
        try:
            cutoff = float(self.box_env_cutoff.text())
            self.state['env_cutoff'] = cutoff
            print(f"Applying envelope cutoff: {cutoff} Hz")
            self._update_plots()
        except ValueError:
            pass  # ignore invalid input
        
    def on_station_changed(self, station_name):
        """
        Handle change of selected station for ZNE view.
        """
        print(f"Changing ZNE station to {station_name}")
        self.state['zne_station'] = station_name
        self._update_plots()

    def _parse_event_name(self, event, idx) -> str:
        return f"Event {idx}: {event['start_time'].strftime('%Y-%m-%d %H:%M')} - {event['end_time'].strftime('%H:%M')}"
    
    def _parse_event_index(self, event_name) -> int:
        return int(event_name.split(":")[0].split(" ")[1])
    
    def on_event_changed(self, event_name):
        # check if self.event_list is empty
        if self.event_list.count() == 0:
            return
        print(f"Changing to event: {event_name}")
        # extract event index from event_name
        event_index = self._parse_event_index(event_name)
        event = self.event_catalog[event_index]
        event_duration = event['end_time'] - event['start_time']
        buffer = 0.05 * event_duration # 5% buffer on each side
        self.state['start_time'] = event['start_time'] - buffer
        self.state['end_time'] = event['end_time'] + buffer
        self.state['current_date'] = UTCDateTime(self.state['start_time'].year,
                                                 self.state['start_time'].month,
                                                 self.state['start_time'].day)
        self.state['selection_start'] = event['start_time']
        self.state['selection_end'] = event['end_time']
        self._update_plots()

    def on_event_added(self):
        # add current selection as new event to catalog
        if self.state['selection_start'] and self.state['selection_end']:
            dialog = EventNoteDialog(self.state['selection_start'], self.state['selection_end'], self)

            if dialog.exec() == QtWidgets.QDialog.Accepted:
                start = dialog.start_time()
                end   = dialog.end_time()
                note  = dialog.note()
                
                new_index = len(self.event_catalog)
                new_event = {
                    'event_id': new_index,
                    'start_time': start,
                    'end_time': end,
                    'notes': note,
                }

                self.event_catalog[len(self.event_catalog)] = new_event
                self.save_temp_catalog()
                #write_event_catalog(self.event_catalog, self.event_catalog_path)
                print(f"Added new event: {new_event}")
                # update event list
                self.event_list.addItem(self._parse_event_name(new_event, new_index))
                # update
    
    def catalog_dialog(self):
        """Ask whether to load or save."""
        diaglog = CatalogDialog(self)
        diaglog.exec()

        if diaglog.choice == CatalogDialog.LOAD:
            self.load_catalog()
        elif diaglog.choice == CatalogDialog.SAVE:
            self.save_catalog()
        elif diaglog.choice == CatalogDialog.NEW:
            confirm = QtWidgets.QMessageBox.question(
                self,
                "Confirm Creating New Catalog",
                "Are you sure you want to create a new catalog? This will erase the current catalog in memory.",
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No
            )
            if confirm == QtWidgets.QMessageBox.Yes:
                self.event_catalog = {}
                self.event_list.clear()
                print("Created new empty catalog.")
        else:
            print("Catalog action cancelled.")

    def load_catalog(self):
        if self.state['last_catalog_load_folder']:
            start_path = os.path.join(self.state['last_catalog_load_folder'], '')
        else:
            start_path = 'data/events'
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Open Catalog",
            start_path,
            "Catalog Files (*.json *.csv *.txt);;All Files (*)"
        )
        if filename:
            print(f"Loading catalog from: {filename}")
            self.event_catalog = read_event_catalog(filename)
            # update event list
            event_names = [
                self._parse_event_name(e, i)
                for i, e in self.event_catalog.items()
            ]
            print(event_names)
            self.event_list.clear()
            self.event_list.addItems(event_names)
            self.state['last_catalog_load_folder'] = os.path.dirname(filename)
            self.state['last_loaded_catalog_file'] = filename
            self.state['all_catalog_changes_saved'] = True
            self._save_state(name='last')

    def save_temp_catalog(self):
        temp_path = os.path.join('data/temp_files', 'temp_catalog.csv')
        write_event_catalog(self.event_catalog, temp_path)
        self.state['all_catalog_changes_saved'] = False
        self._save_state(name='last')

    def save_catalog(self):
        default_name =  self.state.get('last_saved_catalog_filename', f'event_catalog')
        if self.state['last_catalog_save_folder']:
            start_path = os.path.join(self.state['last_catalog_save_folder'], default_name)
        else:
            start_path = os.path.join("outputs", default_name)
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save Catalog",
            start_path,
            "CSV (*.csv);;JSON (*.json);;PyTorch (*.pt *.pth)"
        )
        if filename:
            print(f"Saving catalog to: {filename}")
            write_event_catalog(self.event_catalog, filename)
            self.state['last_catalog_save_folder'] = os.path.dirname(filename)
            self.state['last_saved_catalog_filename'] = os.path.basename(filename)
            self.state['all_catalog_changes_saved'] = True
            self._save_state(name='last')
    
    def _plot_at_scene_pos(self, scene_pos):
        """
        Returns (plot_widget, viewbox, data_point) or (None, None, None)
        """
        for plot in self.stack.currentWidget().plots:
            vb = plot.vb
            if vb.sceneBoundingRect().contains(scene_pos):
                data_pos = vb.mapSceneToView(scene_pos)
                return plot, vb, data_pos
            
        try:
            for plot in self.stack.currentWidget().spectrogram_plots:
                vb = plot.vb
                if vb.sceneBoundingRect().contains(scene_pos):
                    data_pos = vb.mapSceneToView(scene_pos)
                    return plot, vb, data_pos
        except AttributeError:
            pass

        return None, None, None

    def _on_scene_mouse_clicked(self, event):
        # double mouse click event to clear selection
        scene_pos = event.scenePos()
        plot, viewbox, data_pos = self._plot_at_scene_pos(scene_pos)

        if plot is None:
            # Click was NOT on a plot (controls, margins, empty space)
            return

        # Click was on a plot
        if event.double():
            for p in self._get_all_plots():
                vb = p.getViewBox()
                vb.clear_selection()
                vb.clear_permanent_selection()
            self.state['selection_start'] = None
            self.state['selection_end'] = None
        else:
            pass
    
    def edit_catalog(self):
        editor = CatalogEditor(self.event_catalog, parent=self)
        if editor.exec()== QtWidgets.QDialog.Accepted:
            print("Catalog edited.")
            self.event_catalog = editor.get_catalog()
            for k, v in self.event_catalog.items():
                print(f"Event {k}: {v}")
            event_names = [
                self._parse_event_name(e, i)
                for i, e in self.event_catalog.items()
            ]
            self.event_list.clear()
            self.event_list.addItems(event_names)
            self.save_temp_catalog()

if __name__ == "__main__":

    # Argument parsing
    # data directory, stations, fs, max_npts, clip_amplitude, zero_amplitude, init_state
    parser = None

    #data_dir = '/Users/kianhunziker/Documents/UNI/UNIBAS/MA/seisLM/data/tremor/Tremor_daily'
    #data_dir = '/Users/kianhunziker/Documents/UNI/UNIBAS/MA/seisLM/data/tremor/tremor2017'
    #data_dir = '/Users/kianhunziker/Documents/UNI/UNIBAS/MA/seisLM/data/tremor/tremor2020'
    data_dir = '/Users/kianhunziker/Documents/UNI/UNIBAS/MA/seisLM/data/tremor/tremor_temp'
    #data_dir = 'data/iris'
    
    
    app = QtWidgets.QApplication(sys.argv)
    pg.setConfigOptions(antialias=True)
    win = TremorViewer(
        data_dir=data_dir,
        stations=['JUDI', 'LAFE', 'INDI', 'JACO'], #, 'LAFE', 'JACO', 'INDI', ],
        fs=100,
        max_npts=720_000,
        clip_amplitude=20e6,
        zero_amplitude=25e6,
        init_from_state='last',
    )
    win.show()
    sys.exit(app.exec())
