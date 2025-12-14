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


class ThreeComponentSingleStationView(BaseSeismicView):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        self.time_axis = TimeAxis(start_time=self.state['start_time'])
        self.build_plots()


    def build_plots(self):
        self.clear()

        duration = self.state['end_time'] - self.state['start_time']
        self.t = np.arange(int(duration * self.fs + 1)) / self.fs  # Assuming 50 Hz sample rate
        self.env_t = env_t = np.linspace(self.t[0], self.t[-1], len(self.t)//self.fs + 1)

        components = ['Z', 'N', 'E']
        for i, component in enumerate(components):
            vb = YZoomOnlyViewBox()

            if i == len(components) - 1:
                p = pg.PlotItem(
                    axisItems={'bottom': self.time_axis},
                    viewBox=vb
                )
                date_str = self.state['start_time'].strftime("%Y-%m-%d")
                p.setLabel('bottom', f"Station {self.state['zne_station']} on {date_str}")
            else:
                p = pg.PlotItem(viewBox=vb)
            # set label, horizontal rotation
            # TODO: rotate label
            p.setLabel('left', f"{component}")
            p.showGrid(x=True, y=True)

            curve = p.plot([0], [0], pen=pg.mkPen(self.wf_color, width=self.wf_linewidth))
            env = p.plot([0], [0], pen=pg.mkPen(self.env_color, width=self.env_linewidth))

            self.graphics.addItem(p, row=i, col=0)

            self.plots.append(p)
            self.curves.append(curve)
            self.envelopes.append(env)

            if i < len(components) - 1:
                p.getAxis('bottom').setTicks([])  # hide ticks
                p.getAxis('bottom').setStyle(showValues=False)  # hide numbers

        for i in range(1, len(self.plots)):
            self.plots[i].setXLink(self.plots[0])
        # Initial Y-limits
        for p in self.plots:
            p.setYRange(-4, 4)

    
    def refresh(self):
        self.time_axis.start_time = self.state['start_time']

        full_stream = self.stream_dict[self.state['zne_station']]

        for i, component in enumerate(['Z', 'N', 'E']):
            st = select_time_window(
                full_stream,
                self.state['start_time'],
                self.state['end_time']
            )
            st = select_component(st, component)
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
                start_offset = self.state['selection_start'] - self.state['start_time']
                end_offset = self.state['selection_end'] - self.state['start_time']
                vb = self.plots[i].getViewBox()
                if vb.lr_permanent is not None:
                    vb.removeItem(vb.lr_permanent)
                    vb.lr_permanent = None
                vb.lr_permanent = pg.LinearRegionItem(values=(start_offset, end_offset))
                vb.addItem(vb.lr_permanent)
                vb.lr_permanent.setZValue(1000)

        date_str = self.state['start_time'].strftime("%Y-%m-%d")
        self.plots[-1].setLabel('bottom', f"Station {self.state['zne_station']} on {date_str}")
        
        duration = self.t[-1]
        for i, (p, curve, station) in enumerate(zip(self.plots, self.curves, self.stations)):
            # adjust x-range
            duration = len(self.t) / self.fs
            p.setXRange(0, duration)  # convert samples back to seconds


class ThreeComponentMotionPlotView(ThreeComponentSingleStationView):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
    
    def set_square_range(self, plot, x_data, y_data, padding=0.1):
        """
        Adjust plot X and Y range so both axes cover the same range and appear square.
        """
        xmin, xmax = np.min(x_data), np.max(x_data)
        ymin, ymax = np.min(y_data), np.max(y_data)
        
        # determine the total range
        xspan = xmax - xmin
        yspan = ymax - ymin
        span = max(xspan, yspan)  # use the larger span for both axes

        # optional padding
        span *= 1 + padding

        # center
        xmid = 0.5 * (xmax + xmin)
        ymid = 0.5 * (ymax + ymin)

        plot.setXRange(xmid - 0.5*span, xmid + 0.5*span, padding=0)
        plot.setYRange(ymid - 0.5*span, ymid + 0.5*span, padding=0)

    def set_centered_range(self, plot, x_data, y_data, x_padding=0.1, y_padding=0.1):
        """
        Adjust plot X and Y range so both axes are centered around zero with specified padding.
        """
        xmax = np.max(np.abs(x_data))
        ymax = np.max(np.abs(y_data))

        xmax *= 1 + x_padding
        ymax *= 1 + y_padding

        plot.setXRange(-xmax, xmax, padding=0)
        plot.setYRange(-ymax, ymax, padding=0)

    def build_plots(self):
        super().build_plots()

        # add fourth plot 
        # fourth plot is a N-E particle motion plot
        # it is quadratic and shares the x-axis with the other plots
        # it is located on the right side of the other plots
        pm_plot = pg.PlotItem()
        pm_plot.setAspectLocked(True)  # force equal scaling
        pm_plot.showGrid(x=True, y=True)
        pm_curve = pm_plot.plot([0], [0], pen=pg.mkPen(self.wf_color, width=self.wf_linewidth))
        pm_plot.setLabel('bottom', 'E')
        pm_plot.setLabel('left', 'N')

        # Add to the right of waveform plots, spanning all three rows
        self.graphics.addItem(pm_plot, row=0, col=1, rowspan=3)

        # Store particle motion curve for later update
        self.pm_plot = pm_plot
        self.pm_curve = pm_curve
    
    def refresh(self):
        super().refresh()

        # update particle motion plot
        motion_start = self.state['selection_start'] or self.state['start_time']
        motion_end = self.state['selection_end'] or self.state['end_time']
        n_e_data = []
        for c in ['N', 'E']:
            data = self.stream_dict[self.state['zne_station']]
            st = select_time_window(
                data,
                motion_start,
                motion_end
            )
            st = select_component(st, c)
            if self.state['freqmin'] and self.state['freqmax']:
                st = bandpass_filter(
                    st,
                    self.state['freqmin'],
                    self.state['freqmax'],
                    fs=self.fs
                )
            #max_pts = 50_000
            if len(st[0].data) > self.max_npts:
                st[0].data = st[0].data[::len(st[0].data)//self.max_npts]
            n_e_data.append(st[0].data) # no standardization, we want relative amplitudes here

        n_data = n_e_data[0]
        e_data = n_e_data[1]
        self.pm_curve.setData(e_data, n_data)
        # adjust particle motion plot range
        self.set_centered_range(self.pm_plot, n_data, e_data, x_padding=0.1, y_padding=0.1)

        
        