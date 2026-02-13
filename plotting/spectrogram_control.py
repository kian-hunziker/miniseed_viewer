from PySide6 import QtWidgets
from PySide6 import QtCore


class SpectrogramControlWindow(QtWidgets.QWidget):
    paramsChanged = QtCore.Signal(bool)

    def __init__(self, state, parent=None):
        super().__init__(parent)
        self.state = state
        self.setWindowTitle("Spectrogram Settings")
        self.setFixedWidth(300)

        layout = QtWidgets.QFormLayout(self)

        # NFFT spin box in powers of 2
        self.box_nfft = QtWidgets.QSpinBox()
        self.box_nfft.setRange(16, 8192)
        self.box_nfft.setSingleStep(16)
        self.box_nfft.setValue(state['spectrogram']['nperseg'])

        # overlap in the range 0-1
        self.box_overlap = QtWidgets.QDoubleSpinBox()
        self.box_overlap.setRange(0, 1)
        self.box_overlap.setSingleStep(0.1)
        self.box_overlap.setValue(state['spectrogram']['overlap'])

        self.box_fmin = QtWidgets.QDoubleSpinBox()
        self.box_fmin.setRange(0, 100)
        self.box_fmin.setValue(state['spectrogram']['fmin'])

        self.box_fmax = QtWidgets.QDoubleSpinBox()
        self.box_fmax.setRange(0, 100)
        self.box_fmax.setValue(state['spectrogram']['fmax'])

        # cmap selection combo box
        # options: viridis, plasma, inferno, magma, cividis and gray
        self.cmap_box = QtWidgets.QComboBox()
        cmaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'black']
        self.cmap_box.addItems(cmaps)
        current_cmap = state['spectrogram']['cmap']
        if current_cmap in cmaps:
            self.cmap_box.setCurrentText(current_cmap)

        self.independent_filtering_checkbox = QtWidgets.QCheckBox("Independent Bandpass Filtering")
        self.independent_filtering_checkbox.setChecked(state['spectrogram'].get('independent_filtering', False))
        self.independent_filtering_checkbox.stateChanged.connect(self._on_changed)

        self.filter_min_box = QtWidgets.QDoubleSpinBox()
        self.filter_min_box.setRange(0.01, 100)
        self.filter_min_box.setValue(state['spectrogram'].get('f_bandpass_min', 0.1))
        self.filter_min_box.valueChanged.connect(self._on_changed)

        self.filter_max_box = QtWidgets.QDoubleSpinBox()
        self.filter_max_box.setRange(0.1, 100)
        self.filter_max_box.setValue(state['spectrogram'].get('f_bandpass_max', 20.0))
        self.filter_max_box.valueChanged.connect(self._on_changed)

        layout.addRow("NFFT", self.box_nfft)
        layout.addRow("Overlap", self.box_overlap)
        layout.addRow("F min (Hz)", self.box_fmin)
        layout.addRow("F max (Hz)", self.box_fmax)
        layout.addRow("Colormap", self.cmap_box)

        layout.addRow(self.independent_filtering_checkbox)
        layout.addRow("Bandpass F max (Hz)", self.filter_max_box)
        layout.addRow("Bandpass F min (Hz)", self.filter_min_box)
        
        for w in (
            self.box_nfft,
            self.box_overlap,
            self.box_fmin,
            self.box_fmax,
            self.filter_max_box,
            self.filter_min_box,
        ):
            w.valueChanged.connect(self._on_changed)
        # cmap box
        self.cmap_box.currentTextChanged.connect(self._on_changed)
        self.independent_filtering_checkbox.stateChanged.connect(self._on_changed)

    def _on_changed(self):
        spec = self.state['spectrogram']
        spec['nperseg'] = self.box_nfft.value()
        spec['overlap'] = self.box_overlap.value()
        spec['fmin'] = self.box_fmin.value()
        spec['fmax'] = self.box_fmax.value()
        spec['cmap'] = self.cmap_box.currentText()
        spec['independent_filtering'] = self.independent_filtering_checkbox.isChecked()
        self.state['spectrogram']['f_bandpass_min'] = self.filter_min_box.value()
        self.state['spectrogram']['f_bandpass_max'] = self.filter_max_box.value()
        self.paramsChanged.emit(True)
