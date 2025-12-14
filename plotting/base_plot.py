from PySide6 import QtWidgets
from PySide6 import QtCore
import pyqtgraph as pg

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
    