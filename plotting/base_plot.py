from PySide6 import QtWidgets
from PySide6 import QtCore
import pyqtgraph as pg
import numpy as np

from processing.data_processing import (
    bandpass_filter,
    normalize_stream,
    select_time_window,
    select_component,
    compute_envelope,
)

class BaseSeismicView(QtWidgets.QWidget):
    """
    Base class for seismic waveform plotting views.
    Plots in multistation_plots.py and zne_plots.py inherit from this class.
    """
    def __init__(self, 
                 state: dict,
                 stations: list[str],
                 stream_dict: dict,
                 max_npts: int,
                 fs: int,
                 wf_color: str,
                 wf_linewidth: float,
                 env_color: str,
                 env_linewidth: float,
                 ):
        """
        Initialize the base seismic view.
        
        Args:
            state: Application state dict.
            stations: List of station identifiers.
            stream_dict: Dictionary of seismic streams.
            max_npts: Maximum number of points to plot.
            fs: Sampling frequency.
            wf_color: Waveform color.
            wf_linewidth: Waveform line width.
            env_color: Envelope color.
            env_linewidth: Envelope line width.
        """
        super().__init__()
        self.state = state
        self.stream_dict = stream_dict
        self.stations = stations
        self.max_npts = max_npts
        self.fs = fs
        self.wf_color = wf_color
        self.wf_linewidth = wf_linewidth
        self.env_color = env_color
        self.env_linewidth = env_linewidth

        layout = QtWidgets.QVBoxLayout(self)
        self.graphics = pg.GraphicsLayoutWidget()
        layout.addWidget(self.graphics)

        self.plots = []
        self.curves = []
        self.envelopes = []

        self.t = None  # time array for waveforms
        self.env_t = None  # time array for envelopes

        self.shift_pressed = False
        self._install_shift_events()
    
    def _install_shift_events(self):
        # Catch key events from the graphics widget
        self.graphics.keyPressEvent = self._key_press
        self.graphics.keyReleaseEvent = self._key_release
    
    def _key_press(self, event):
        if event.key() == QtCore.Qt.Key_Shift:
            self.shift_pressed = True
            # optionally change cursor
            self.graphics.setCursor(QtCore.Qt.CrossCursor)
        else:
            # call default handler
            super().keyPressEvent(event)

    def _key_release(self, event):
        if event.key() == QtCore.Qt.Key_Shift:
            self.shift_pressed = False
            self.graphics.setCursor(QtCore.Qt.ArrowCursor)

    def clear(self):
        self.graphics.clear()
        self.plots.clear()
        self.curves.clear()
        self.envelopes.clear()

    def refresh(self):
        raise NotImplementedError
    
    def update_wf_plot(self, i, station, component):
        full_stream = self.stream_dict[station]

        st = select_time_window(
            full_stream,
            self.state['start_time'],
            self.state['end_time']
        )
        st = select_component(st, component=component)

        if self.state['freqmin'] and self.state['freqmax']:
            st = bandpass_filter(
                st,
                self.state['freqmin'],
                self.state['freqmax'],
                fs=self.fs
            )

        data = st[0].data
        sample_interval = max(1, len(data) // self.max_npts)
        self.state['downsampling_factor'] = sample_interval

        if i == 0:
            self.t = np.arange(len(data)) / self.fs
            self.t = self.t[::sample_interval] / sample_interval
            self.time_axis.set_resampling_factor(sample_interval)

        if len(data) > self.max_npts:
            data = data[::sample_interval]

        if len(data) > len(self.t):
            data = data[:len(self.t)]
        elif len(data) < len(self.t):
            data = np.pad(data, (0, len(self.t) - len(data)), 'constant')

        data = normalize_stream(data)

        if self.state['show_waveform']:
            self.curves[i].setData(self.t, data)
            # adjust y-range according to std
            self.plots[i].setYRange(-5, 5)
        else:
            self.curves[i].clear()

        if self.state['show_envelope']:
            env = compute_envelope(
                data,
                cutoff=self.state['env_cutoff'],
                fs=self.fs / sample_interval
            )
            env_t = np.linspace(self.t[0], self.t[-1], len(env))
            self.envelopes[i].setData(env_t, env)
        else:
            self.envelopes[i].clear()

        # mark selection region if any
        if self.state['selection_start'] and self.state['selection_end']:
            start_offset = (self.state['selection_start'] - self.state['start_time']) / self.state['downsampling_factor']
            end_offset = (self.state['selection_end'] - self.state['start_time']) / self.state['downsampling_factor']
            vb = self.plots[i].getViewBox()
            if vb.lr_permanent is not None:
                vb.removeItem(vb.lr_permanent)
                vb.lr_permanent = None
            vb.lr_permanent = pg.LinearRegionItem(values=(start_offset, end_offset))
            vb.addItem(vb.lr_permanent)
            vb.lr_permanent.setZValue(1000)

    def get_time_axis_format(self) -> str:
        duration_seconds = self.state['end_time'] - self.state['start_time']
        if duration_seconds <= 600:
            return "%H:%M:%S"
        elif self.state['start_time'].day != self.state['end_time'].day:
            return "%Y-%m-%d\n%H:%M"
        else:
            return "%H:%M"
    