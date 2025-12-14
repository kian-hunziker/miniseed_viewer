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

from processing.data_loading import get_all_miniseed_files, load_waveform_multiple_days

from plotting.mutlistation_plots import MultiStationSingleComponentView
from plotting.zne_plots import ThreeComponentSingleStationView, ThreeComponentMotionPlotView

from processing.data_processing import (
    clip_trace_amplitude,
    set_high_amplitude_gaps_to_zero,
)


# ==============================================================================
# MAIN APPLICATION
# ==============================================================================

class TremorViewer(QtWidgets.QWidget):
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
        
        # load waveform data for all stations
        self.waveform_files = get_all_miniseed_files(data_dir=self.data_dir)
        self.station_streams = {}
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
            st = set_high_amplitude_gaps_to_zero(st, amplitude_threshold=zero_amplitude)
            st = clip_trace_amplitude(st, max_amplitude=clip_amplitude)
            self.station_streams[station] = st
        
        # setup window
        self.setWindowTitle("Tremor Viewer v3.0")
        self.resize(1200, 600)

        # initial state
        if init_from_state in ['default', 'last']:
            try:
                self._load_state(name=init_from_state)
            except:
                print(f"Failed to load state from config/{init_from_state}.json, using default state.")
                self._setup_initial_state()
        else:
            self._setup_initial_state()
            
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

        layout = QtWidgets.QVBoxLayout(self)

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
        self.all_plots = self.multi_station_view.plots + self.three_comp_view.plots + self.three_comp_with_motion_view.plots

        # default view
        if self.state['current_view'] == 'three_comp' and not self.state['show_motion_plot']:
            self.stack.setCurrentWidget(self.three_comp_view)
        elif self.state['current_view'] == 'three_comp' and self.state['show_motion_plot']:
            self.stack.setCurrentWidget(self.three_comp_with_motion_view)
        elif self.state['current_view'] == 'multi_station':
            self.stack.setCurrentWidget(self.multi_station_view)
        else:
            raise ValueError(f"Invalid current_view: {self.state['current_view']}")
        layout.addWidget(self.stack)
        
        # build controls
        self._build_controls(layout)
        
        # initial refresh
        QtCore.QTimer.singleShot(0, self.stack.currentWidget().refresh)
        
        # setup dragging signals
        self._setup_dragging_signals()

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
                'last_save_folder': None,
            }
        # override saved default state, to update defaults
        self._save_state(name='default')
        
    
    def _get_all_plots(self):
        return self.all_plots
    
    def _setup_dragging_signals(self):
        for p in self._get_all_plots():
            vb = p.getViewBox()
            vb.dragging_signal.connect(self.on_dragging_changed)
    
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
        
        if self.state['current_view'] == 'multi_station':
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

        # envelope and waveform visibility
        self.cb_show_waveform = QtWidgets.QCheckBox("Show Waveform")
        self.cb_show_waveform.setChecked(self.state['show_waveform'])
        controls.addWidget(self.cb_show_waveform, 0, 5)
        self.cb_show_envelope = QtWidgets.QCheckBox("Show Envelope")
        self.cb_show_envelope.setChecked(self.state['show_envelope'])
        controls.addWidget(self.cb_show_envelope, 0, 6)
        self.cb_show_motion_plot = QtWidgets.QCheckBox("Show Motion Plot")
        self.cb_show_motion_plot.setChecked(self.state['show_motion_plot'])
        controls.addWidget(self.cb_show_motion_plot, 0, 7)
        
        self.box_env_cutoff = QtWidgets.QLineEdit(str(self.state['env_cutoff']))
        controls.addWidget(self.box_env_cutoff, 1, 5)
        self.btn_apply_env_cutoff = QtWidgets.QPushButton("Apply Env Cutoff")
        controls.addWidget(self.btn_apply_env_cutoff, 1, 6)

        # Reset Y button
        self.btn_reset_y = QtWidgets.QPushButton("Reset Y-Limits")
        controls.addWidget(self.btn_reset_y, 4, 4)

        # save button
        self.btn_save = QtWidgets.QPushButton("Save Figure")
        controls.addWidget(self.btn_save, 4, 7)

        # Connect signals
        self.station_selector.currentTextChanged.connect(self.on_station_changed)
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
        current_view = self.stack.currentWidget()
        current_view.refresh()

    def on_dragging_changed(self, dragging):
        # update dragging state
        self.state['dragging'] = dragging

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
            else:
                self.stack.setCurrentWidget(self.three_comp_view)
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
        if self.state['last_save_folder']:
            start_path = os.path.join(self.state['last_save_folder'], default_name)
        else:
            start_path = os.path.join("outputs", default_name)
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            None, "Save plot", start_path, "PNG (*.png);;SVG (*.svg)"
        )
        # update last save folder
        if path:
            self.state['last_save_folder'] = os.path.dirname(path)
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
    

if __name__ == "__main__":

    # Argument parsing
    # data directory, stations, fs, max_npts, clip_amplitude, zero_amplitude, init_state
    parser = None

    data_dir = '/Users/kianhunziker/Documents/UNI/UNIBAS/MA/seisLM/data/tremor/Tremor_daily'
    
    
    app = QtWidgets.QApplication(sys.argv)
    pg.setConfigOptions(antialias=True)
    win = TremorViewer(
        data_dir=data_dir,
        stations=['JUDI', 'LAFE', 'JACO', 'INDI', ],
        fs=100,
        max_npts=720_000,
        clip_amplitude=20e6,
        zero_amplitude=25e6,
        init_from_state='last',
    )
    win.show()
    sys.exit(app.exec())
