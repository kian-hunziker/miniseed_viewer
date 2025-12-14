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
            #data = np.random.randn(len(self.t))  # placeholder
            #print(data)

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

        for i, station in enumerate(self.stations):
            full_stream = self.stream_dict[station]

            st = select_time_window(
                full_stream,
                self.state['start_time'],
                self.state['end_time']
            )
            st = select_component(st, self.state['component'])

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

            data = normalize_stream(data)

            if self.state['show_waveform']:
                self.curves[i].setData(self.t, data)
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
        
        date_str = self.state['start_time'].strftime("%Y-%m-%d")
        self.plots[-1].setLabel('bottom', date_str)

        duration = self.t[-1]
        for i, (p, curve, station) in enumerate(zip(self.plots, self.curves, self.stations)):
            # adjust x-range
            duration = len(self.t) / self.fs
            p.setXRange(0, duration)  # convert samples back to seconds
            #p.setYRange(self.state['y_limits'][station][0], self.state['y_limits'][station][1])

