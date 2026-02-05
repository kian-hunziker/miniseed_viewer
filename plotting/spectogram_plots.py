import numpy as np
import pyqtgraph as pg
import numpy as np

from plotting.base_plot import BaseSeismicView
from plotting.mutlistation_plots import MultiStationSingleComponentView
from plotting.utility_classes import TimeAxis, YZoomOnlyViewBox

from processing.data_processing import (
    bandpass_filter,
    normalize_stream,
    select_time_window,
    select_component,
    compute_envelope,
)
from scipy.signal import spectrogram


class MultistationWithSpectogram(BaseSeismicView):
    def __init__(self, state, **kwargs):
        super().__init__(state, **kwargs)
        self.time_axis = TimeAxis(start_time=self.state['start_time'])
        self.time_axis_left = TimeAxis(start_time=self.state['start_time'])
        self.spectrogram_plots = []
        self.spectrogram_images = []
        self.color_map = 'viridis'
        self.build_plots()

    def clear(self):
        self.graphics.clear()
        self.plots.clear()
        self.curves.clear()
        self.envelopes.clear()
        self.spectrogram_plots.clear()
        self.spectrogram_images.clear()

    def update_color_map(self, cmap_name):
        self.color_map = cmap_name
        self.spec_map = pg.colormap.get(cmap_name)
        for img in self.spectrogram_images:
            img.setLookupTable(self.spec_map.getLookupTable(0.0, 1.0, 256))

    def build_plots(self):
        self.clear()

        # build waveform plots
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
            

    
        # Add spectrograms below each station plot
        self.spec_map = pg.colormap.get(self.color_map)
        if self.state['spectrogram']['cmap'] != self.color_map:
            self.update_color_map(self.state['spectrogram']['cmap'])
        for i, station in enumerate(self.stations):
            vb = YZoomOnlyViewBox()

            if i == len(self.stations) - 1:
                p = pg.PlotItem(
                    axisItems={'bottom': self.time_axis_left},
                    viewBox=vb
                )
            else:
                p = pg.PlotItem(viewBox=vb)

            p.setLabel('left', f"{station}")
            p.showGrid(x=True, y=True)

            spectrogram = pg.ImageItem()
            spectrogram.setLookupTable(self.spec_map.getLookupTable(0.0, 1.0, 256))
            p.addItem(spectrogram)

            self.graphics.addItem(p, row=i, col=1)

            # Link x-axes
            p.setXLink(self.plots[i])

            self.spectrogram_plots.append(p)
            self.spectrogram_images.append(spectrogram)

            if i < len(self.stations) - 1:
                p.getAxis('bottom').setStyle(showValues=False)


    def update_spectrogram(self, i, station):
        full_stream = self.stream_dict[station]

        st = select_time_window(
            full_stream,
            self.state['start_time'],
            self.state['end_time']
        )
        st = select_component(st, component=self.state['component'])

        if self.state['freqmin'] and self.state['freqmax']:
            st = bandpass_filter(
                st,
                self.state['freqmin'],
                self.state['freqmax'],
                fs=self.fs
            )

        data = st[0].data
        data = normalize_stream(data)

        spectrogram_params = self.state['spectrogram']
        mode = spectrogram_params['mode']
        nperseg = spectrogram_params['nperseg']
        overlap = spectrogram_params['overlap']
        noverlap = int(nperseg * overlap)
        print(f'noverlap: {noverlap}')
        
        f, t_spec, Sxx = spectrogram(
            data, 
            fs=self.fs, 
            nperseg=nperseg, 
            noverlap=noverlap,
            mode=mode
            )

        if spectrogram_params['log_scale']:
            Sxx = np.log10(Sxx + 1e-12)
        Sxx = Sxx.T
        # Update spectrogram image
        img = self.spectrogram_images[i]

        #vmin, vmax = np.percentile(Sxx, [5, 95])
        #img.setLevels([vmin, vmax])
        Sxx -= Sxx.min()
        Sxx /= Sxx.max()
        img.setImage(Sxx, autoLevels=True) # if set to false, it seems to crash

        dt = t_spec[1] - t_spec[0]
        df = f[1] - f[0]

        tr = pg.QtGui.QTransform()
        tr.scale(dt, df)

        img.setTransform(tr)
        img.setPos(t_spec[0], f[0])

        #self.spectrogram_images[i].scale(t_spec[1] - t_spec[0], f[1] - f[0])
        self.spectrogram_images[i].setPos(t_spec[0], f[0])
        #self.spectrogram_plots[i].setLabel('bottom', 'Time (s)')
        #self.spectrogram_plots[i].setLabel('left', 'Frequency (Hz)')

        # mark selection region if any
        if self.state['selection_start'] and self.state['selection_end']:
            start_offset = (self.state['selection_start'] - self.state['start_time']) / self.state['downsampling_factor']
            end_offset = (self.state['selection_end'] - self.state['start_time']) / self.state['downsampling_factor']
            vb = self.spectrogram_plots[i].getViewBox()
            if vb.lr_permanent is not None:
                vb.removeItem(vb.lr_permanent)
                vb.lr_permanent = None
            vb.lr_permanent = pg.LinearRegionItem(values=(start_offset, end_offset))
            vb.addItem(vb.lr_permanent)
            vb.lr_permanent.setZValue(1000)
        
    def refresh(self):
        # Refresh time axis
        self.time_axis.start_time = self.state['start_time']
        self.time_axis.set_format(self.get_time_axis_format())
        self.time_axis_left.start_time = self.state['start_time']
        self.time_axis_left.set_format(self.get_time_axis_format())

        for i, station in enumerate(self.stations):
            self.update_wf_plot(i, station, self.state['component'])
            self.update_spectrogram(i, station)
        date_str = self.state['start_time'].strftime("%Y-%m-%d")
        self.plots[-1].setLabel('bottom', date_str)
        self.spectrogram_plots[-1].setLabel('bottom', date_str)

        duration = len(self.t) / self.fs
        for p in self.plots:
            # adjust x-range
            # this is corrected for downsampling in TimeAxis
            p.setXRange(0, duration)  # convert samples back to seconds
        print("Updated spectrogram plots")

