#!/usr/bin/env python3
import sys
import numpy as np
from PySide6 import QtWidgets
import pyqtgraph as pg
import os
import glob
import numpy as np
from scipy.signal import hilbert
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons
from matplotlib.widgets import TextBox, Button
from obspy import UTCDateTime, Stream, Trace
from pyqtgraph.exporters import ImageExporter, SVGExporter
from tqdm import tqdm

from processing.data_loading import get_all_miniseed_files, load_waveform_multiple_days

from datetime import datetime, timedelta
from processing.data_processing import (
    clip_trace_amplitude,
    set_high_amplitude_gaps_to_zero,
    bandpass_filter,
    normalize_stream,
    select_time_window,
    select_component,
    compute_envelope,
)



# ==============================================================================
# Prepare Data
# ==============================================================================
FS = 100  # sample rate in Hz
MAX_NPTS = 720_000

waveform_files = get_all_miniseed_files()
stations = ['JUDI', 'LAFE', 'JACO', 'INDI', ]#'JACO',
station_streams = {}
for station in tqdm(stations, desc="Loading stations"):
    st = load_waveform_multiple_days(
        station=station,
        waveform_files=waveform_files,
        fs=FS,
    )
    station_streams[station] = st

clip_threshold = 1e6
zero_threshold = 25e6
for station in stations:
    st = station_streams[station]
    st = set_high_amplitude_gaps_to_zero(st, amplitude_threshold=zero_threshold)
    st = clip_trace_amplitude(st, max_amplitude=clip_threshold)
    station_streams[station] = st

# ==============================================================================
# Auxiliary Plot Classes
# ==============================================================================

class TimeAxis(pg.AxisItem):
    def __init__(self, start_time, **kwargs):
        super().__init__(orientation='bottom', **kwargs)
        self.start_time = start_time  # should be a datetime object
        self.resampling_factor = 1

    def tickStrings(self, values, scale, spacing):
        strings = []
        for x in values:
            t = self.start_time + self.resampling_factor * timedelta(seconds=x)
            strings.append(t.strftime("%H:%M"))
        return strings
    
    def set_resampling_factor(self, factor):
        self.resampling_factor = factor
    
    
    
class YZoomOnlyViewBox(pg.ViewBox):
    """Custom ViewBox that only zooms in Y, leaves X unchanged."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMouseMode(self.PanMode)  # horizontal pan works normally

    def wheelEvent(self, ev, axis=None):
        """Zoom only Y, X remains fixed."""
        factor = 0.9 if ev.delta() > 0 else 1.1
        y_min, y_max = self.viewRange()[1]
        center = 0.5 * (y_min + y_max)
        half = 0.5 * (y_max - y_min) * factor
        self.setYRange(center - half, center + half, padding=0)
        ev.accept()  # stop propagation, prevents default X+Y zoom


# ==============================================================================
# MAIN APPLICATION
# ==============================================================================

class IndependentYZoomPlots(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Tremor Viewer v3.0")
        self.resize(1200, 600)
        self.stations = stations
        self.state = {
            'component': 'N',
            'start_time': UTCDateTime(2020, 8, 24, 5, 0, 0),
            'end_time': UTCDateTime(2020, 8, 24, 7, 0, 0),
            'current_date': UTCDateTime(2020, 8, 24),
            'freqmin': None,
            'freqmax': None,
            'y_limits': {station: (-4, 4) for station in stations},
            'show_envelope': False,
            'show_waveform': True,
            'env_cutoff': 0.1,
        }
        self.wf_color = 'k'#(100, 100, 100)  # waveform color
        self.wf_linewidth = 0.8  # waveform line width
        self.env_color = 'r'  # envelope color
        self.env_linewidth = 2  # envelope line width
        self.background_color = 'w' #(180, 180, 180)  # dark background

        self.last_save_folder = None

        layout = QtWidgets.QVBoxLayout(self)
        pg.setConfigOption('background', self.background_color)
        pg.setConfigOption('foreground', self.wf_color)
        pg.setConfigOptions(useOpenGL=True)
        self.graphics = pg.GraphicsLayoutWidget()
        layout.addWidget(self.graphics)

        # Controls: grid
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
        self.box_date = QtWidgets.QLineEdit("08-24")
        controls.addWidget(self.box_date, 1, 1)
        self.btn_apply_date = QtWidgets.QPushButton("Apply Date")
        controls.addWidget(self.btn_apply_date, 1, 2)

        # Bandpass controls
        controls.addWidget(QtWidgets.QLabel("Fmin (Hz):"), 2, 0)
        self.box_fmin = QtWidgets.QLineEdit("2")
        controls.addWidget(self.box_fmin, 2, 1)
        controls.addWidget(QtWidgets.QLabel("Fmax (Hz):"), 2, 2)
        self.box_fmax = QtWidgets.QLineEdit("8")
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
        self.rbN.setChecked(True)
        comp_layout.addWidget(self.rbZ)
        comp_layout.addWidget(self.rbN)
        comp_layout.addWidget(self.rbE)

        # envelope and waveform visibility
        self.cb_show_waveform = QtWidgets.QCheckBox("Show Waveform")
        self.cb_show_waveform.setChecked(self.state['show_waveform'])
        controls.addWidget(self.cb_show_waveform, 0, 5)
        self.cb_show_envelope = QtWidgets.QCheckBox("Show Envelope")
        self.cb_show_envelope.setChecked(self.state['show_envelope'])
        controls.addWidget(self.cb_show_envelope, 0, 6)
        
        self.box_env_cutoff = QtWidgets.QLineEdit(str(self.state['env_cutoff']))
        controls.addWidget(self.box_env_cutoff, 1, 5)
        self.btn_apply_env_cutoff = QtWidgets.QPushButton("Apply Env Cutoff")
        controls.addWidget(self.btn_apply_env_cutoff, 1, 6)

        # Reset Y button
        self.btn_reset_y = QtWidgets.QPushButton("Reset Y-Limits")
        controls.addWidget(self.btn_reset_y, 4, 4)

        # save button
        self.btn_save = QtWidgets.QPushButton("Save Figure")
        controls.addWidget(self.btn_save, 4, 3)

        # time axis
        # in seconds after start_time
        duration = self.state['end_time'] - self.state['start_time']
        self.t = np.arange(int(duration * FS + 1)) / FS  # Assuming 50 Hz sample rate
        self.env_t = env_t = np.linspace(self.t[0], self.t[-1], len(self.t)//FS + 1)

        # PLOTS FOR WAVEFORMS
        self.plots = []
        self.curves = []
        self.envelopes = []
        self.time_axis = TimeAxis(start_time=self.state['start_time'])
        for i, station in enumerate(stations):
            # Prepare data
            data = station_streams[station]
            data = select_time_window(data, self.state['start_time'], self.state['end_time'])
            data = select_component(data, self.state['component'])
            data = normalize_stream(data[0].data)

            # Create plot with custom ViewBox
            vb1 = YZoomOnlyViewBox()
            if i == len(stations) - 1:
                p = pg.PlotItem(axisItems={'bottom': self.time_axis}, viewBox=YZoomOnlyViewBox())
                date_str = self.state['start_time'].strftime("%Y-%m-%d")
                p.setLabel('bottom', date_str)
            else:
                p = pg.PlotItem(viewBox=vb1)
            
            # plot waveform and envelope
            curve = p.plot(self.t, data, pen=pg.mkPen(self.wf_color, width=self.wf_linewidth))
            if self.state['show_envelope']:
                env = p.plot(self.t[::FS], compute_envelope(data), pen=pg.mkPen(self.env_color, width=self.env_linewidth))  # red envelope
            else:
                env = p.plot([], [], pen=pg.mkPen(self.env_color, width=self.env_linewidth))  # empty envelope if not shown
            p.setLabel('left', station)
            p.showGrid(x=True, y=True)
            
            self.graphics.addItem(p, row=i, col=0)
            self.plots.append(p)
            self.curves.append(curve)
            self.envelopes.append(env)
            # Only show x-axis labels on the last plot
            if i < len(stations) - 1:
                p.getAxis('bottom').setTicks([])  # hide ticks
                p.getAxis('bottom').setStyle(showValues=False)  # hide numbers


        # Link all X axes (optional: horizontal panning)
        for i in range(1, len(self.plots)):
            self.plots[i].setXLink(self.plots[0])

        # Initial Y-limits
        for p, station in zip(self.plots, stations):
            p.setYRange(self.state['y_limits'][station][0], self.state['y_limits'][station][1])

        self.format_x_labels()

        # Connect signals
        self.btn_apply_bp.clicked.connect(self.on_apply_bandpass)
        self.btn_apply_time.clicked.connect(self.on_apply_time)
        self.btn_apply_date.clicked.connect(self.on_apply_date)
        self.rbZ.toggled.connect(self.on_component_changed)
        self.rbN.toggled.connect(self.on_component_changed)
        self.rbE.toggled.connect(self.on_component_changed)
        self.btn_reset_y.clicked.connect(self.on_reset_y)
        self.btn_save.clicked.connect(self.save_plot)
        self.cb_show_envelope.stateChanged.connect(self.show_envelope_changed)
        self.cb_show_waveform.stateChanged.connect(self.show_waveform_changed)
        self.btn_apply_env_cutoff.clicked.connect(self.on_apply_env_cutoff)

        # connect return pressed in text boxes to apply functions
        self.box_start.returnPressed.connect(self.on_apply_time)
        self.box_end.returnPressed.connect(self.on_apply_time)
        self.box_date.returnPressed.connect(self.on_apply_date)
        self.box_fmin.returnPressed.connect(self.on_apply_bandpass)
        self.box_fmax.returnPressed.connect(self.on_apply_bandpass)
        self.box_env_cutoff.returnPressed.connect(self.on_apply_env_cutoff)

    def format_x_labels(self):
        # format x-axis labels as HH:MM based on start_time
        p = self.plots[-1]
        start_time = self.state['start_time']

    def on_reset_y(self):
        # loop through plots and reset y-limits
        for p, station in zip(self.plots, self.stations):
            p.setYRange(self.state['y_limits'][station][0], self.state['y_limits'][station][1])

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
    
    def on_apply_time(self):
        t0 = self._parse_time_input(self.box_start.text())
        t1 = self._parse_time_input(self.box_end.text())
        if t0 is None or t1 is None:
            print("Could not parse times. Use ISO or HH:MM.")
            return
        if t1 <= t0:
            print("End must be after start")
            return
        self.state["current_date"] = UTCDateTime(t0.year, t0.month, t0.day)
        self.state["start_time"] = t0
        self.state["end_time"] = t1
        self.update_plots()

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
        self.update_plots()

    def update_plots(self):
        # loop through stations and update plots
        self.time_axis.start_time = self.state['start_time']
        for i, station in enumerate(self.stations):
            full_stream = station_streams[station]
            # trim to start/end time
            trimmed_stream = select_time_window(full_stream, self.state['start_time'], self.state['end_time'])
            # select component
            selected_stream = select_component(trimmed_stream, self.state['component'])
            # apply bandpass if set
            if self.state['freqmin'] is not None and self.state['freqmax'] is not None:
                filtered_stream = bandpass_filter(selected_stream, self.state['freqmin'], self.state['freqmax'], fs=FS)
            else:
                filtered_stream = selected_stream
            
            # resample if too long
            data = filtered_stream[0].data
            sample_interval = len(data) // MAX_NPTS
            if i == 0:
                self.t = np.arange(len(filtered_stream[0].data)) / FS  # update time axis based on actual data length
                self.t = self.t[::sample_interval] / sample_interval  # adjust time axis accordingly
                self.time_axis.set_resampling_factor(sample_interval)
            if len(data) > MAX_NPTS:
                # decimate data to fit
                data = data[::sample_interval]
                if i == 0:
                    print(f"Data too long ({len(data)} points), truncating to {MAX_NPTS} points.")
            
            # normalize            
            data = normalize_stream(data)

            if self.state['show_waveform']:
                self.curves[i].setData(self.t, data)
            else:
                self.curves[i].setData([], [])  # clear waveform if not shown
            if self.state['show_envelope']:
                env_data = compute_envelope(data, cutoff=self.state['env_cutoff'], fs=FS / sample_interval)
                if i == 0:
                    self.env_t = np.linspace(self.t[0], self.t[-1], len(env_data))
                self.envelopes[i].setData(self.env_t, env_data)
            else:
                self.envelopes[i].setData([], [])  # clear envelope if not shown
            date_str = self.state['start_time'].strftime("%Y-%m-%d")
            self.plots[-1].setLabel('bottom', date_str)
        
        for i, (p, curve, station) in enumerate(zip(self.plots, self.curves, self.stations)):
            # adjust x-range
            duration = len(self.t) / FS
            p.setXRange(0, duration)  # convert samples back to seconds
    
    def save_plot(self):
        default_name = f'tremor_plot_{self.state["current_date"].strftime("%Y%m%d")}_{self.state["component"]}_{self.state["freqmin"]}-{self.state["freqmax"]}Hz'
        if self.last_save_folder:
            start_path = os.path.join(self.last_save_folder, default_name)
        else:
            start_path = os.path.join("notes/tremor_plot", default_name)
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            None, "Save plot", start_path, "PNG (*.png);;SVG (*.svg)"
        )
        # update last save folder
        if path:
            self.last_save_folder = os.path.dirname(path)
        if not path:
            return

        if path.endswith(".png"):
            ImageExporter(self.graphics.scene()).export(path)
            print(f"Saving plot to {path}.png")
        elif path.endswith(".svg"):
            SVGExporter(self.graphics.scene()).export(path)
            print(f"Saving plot to {path}.svg")

    def on_apply_bandpass(self):
        try:
            fmin = float(self.box_fmin.text()) if self.box_fmin.text().strip() != "" else None
            fmax = float(self.box_fmax.text()) if self.box_fmax.text().strip() != "" else None
            self.state['freqmin'] = fmin
            self.state['freqmax'] = fmax
            print(f"Applying bandpass filter: {fmin} - {fmax} Hz")
            self.update_plots()
        except ValueError:
            pass  # ignore invalid input
    
    def on_component_changed(self):
        if self.rbZ.isChecked():
            comp = "Z"
        elif self.rbN.isChecked():
            comp = "N"
        else:
            comp = "E"
        # only change and redraw if actually different
        if comp != self.state["component"]:
            print(f"Changing component to {comp}")
            self.state["component"] = comp
            self.update_plots()

    # mouse release event to print current x-limits
    def mouseReleaseEvent(self, event):
        # get current x-limits from first plot
        x_min, x_max = self.plots[0].viewRange()[0]
        start_time = self.state['start_time'] + x_min
        end_time = self.state['start_time'] + x_max
        self.state['start_time'] = start_time
        self.state['end_time'] = end_time
        self.state['current_date'] = UTCDateTime(start_time.year, start_time.month, start_time.day)
        self.update_plots()
        super().mouseReleaseEvent(event)

    def show_envelope_changed(self, state):
        """Toggle envelope visibility."""
        print(f"Show envelope: {state > 0}")
        self.state['show_envelope'] = state > 0
        self.update_plots()
    
    def on_apply_env_cutoff(self):
        try:
            cutoff = float(self.box_env_cutoff.text())
            self.state['env_cutoff'] = cutoff
            print(f"Applying envelope cutoff: {cutoff} Hz")
            self.update_plots()
        except ValueError:
            pass  # ignore invalid input
    
    def show_waveform_changed(self, state):
        """Toggle waveform visibility."""
        print(f"Show waveform: {state > 0}")
        self.state['show_waveform'] = state > 0
        self.update_plots()

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    pg.setConfigOptions(antialias=True)
    win = IndependentYZoomPlots()
    win.show()
    sys.exit(app.exec())
