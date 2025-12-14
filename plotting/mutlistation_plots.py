import numpy as np
import pyqtgraph as pg
import numpy as np

from plotting.base_plot import BaseSeismicView
from plotting.utility_classes import TimeAxis, YZoomOnlyViewBox

from processing.data_processing import (
    bandpass_filter,
    normalize_stream,
    select_time_window,
    select_component,
    compute_envelope,
)


class MultiStationSingleComponentView(BaseSeismicView):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        self.time_axis = TimeAxis(start_time=self.state['start_time'])
        self.build_plots()

    def build_plots(self):
        self.clear()

        duration = self.state['end_time'] - self.state['start_time']
        self.t = np.arange(int(duration * self.fs + 1)) / self.fs  # Assuming 50 Hz sample rate
        self.env_t = env_t = np.linspace(self.t[0], self.t[-1], len(self.t)//self.fs + 1)

        for i, station in enumerate(self.stations):
            vb = YZoomOnlyViewBox()

            if i == len(self.stations) - 1:
                p = pg.PlotItem(
                    axisItems={'bottom': self.time_axis},
                    viewBox=vb
                )
            else:
                p = pg.PlotItem(viewBox=vb)

            p.setLabel('left', station)
            p.showGrid(x=True, y=True)

            curve = p.plot([0], [0], pen=pg.mkPen(self.wf_color, width=self.wf_linewidth))
            env = p.plot([0], [0], pen=pg.mkPen(self.env_color, width=self.env_linewidth))

            self.graphics.addItem(p, row=i, col=0)

            self.plots.append(p)
            self.curves.append(curve)
            self.envelopes.append(env)

            if i < len(self.stations) - 1:
                p.getAxis('bottom').setStyle(showValues=False)

        for i in range(1, len(self.plots)):
            self.plots[i].setXLink(self.plots[0])
        # Initial Y-limits
        for p, station in zip(self.plots, self.stations):
            p.setYRange(self.state['y_limits'][f"{station}_{self.state['component']}"][0],
                        self.state['y_limits'][f"{station}_{self.state['component']}"][1])
    
    def refresh(self):
        self.time_axis.start_time = self.state['start_time']
        self.time_axis.set_format(self.get_time_axis_format())

        for i, station in enumerate(self.stations):
            self.update_wf_plot(i, station, self.state['component'])
        
        date_str = self.state['start_time'].strftime("%Y-%m-%d")
        self.plots[-1].setLabel('bottom', date_str)

        duration = len(self.t) / self.fs
        for p in self.plots:
            # adjust x-range
            # this is corrected for downsampling in TimeAxis
            p.setXRange(0, duration)  # convert samples back to seconds
