from PySide6 import QtWidgets
from PySide6 import QtCore
import pyqtgraph as pg
from datetime import datetime, timedelta


# ==============================================================================
# Auxiliary Plot Classes
# ==============================================================================

class TimeAxis(pg.AxisItem):
    def __init__(self, start_time, **kwargs):
        super().__init__(orientation='bottom', **kwargs)
        self.start_time = start_time  # should be a datetime object
        self.resampling_factor = 1
        self.format = "%H:%M"

    def tickStrings(self, values, scale, spacing):
        strings = []
        for x in values:
            t = self.start_time + self.resampling_factor * timedelta(seconds=x)
            strings.append(t.strftime(self.format))
        return strings
    
    def set_resampling_factor(self, factor):
        self.resampling_factor = factor

    def set_format(self, fmt):
        """
        Set the datetime format for the axis labels.
        Choose from: ["%H:%M", "%Y-%m-%d %H:%M", "%Y-%m-%d"]
        """
        self.format = fmt

    
class YZoomOnlyViewBox(pg.ViewBox):
    dragging_signal = QtCore.Signal(bool)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMouseMode(self.PanMode)
        self.lr = None
        self.lr_permanent = None
        self.shift_dragging = False
        self.dragging = False

    def mouseDragEvent(self, ev, axis=None):
        # Shift + drag → selection creation
        if ev.modifiers() == QtCore.Qt.ShiftModifier and ev.button() == QtCore.Qt.LeftButton:
            if self.lr is None:
                self.lr = pg.LinearRegionItem()
                self.addItem(self.lr)
                self.lr.setZValue(1000)

            if not self.shift_dragging:
                pos = self.mapSceneToView(ev.scenePos()).x()
                self.lr.setRegion((pos, pos))
                self.shift_dragging = True

            pos = self.mapSceneToView(ev.scenePos()).x()
            r0, r1 = self.lr.getRegion()
            if abs(pos - r0) < abs(pos - r1):
                self.lr.setRegion((pos, r1))
            else:
                self.lr.setRegion((r0, pos))

            ev.accept()
            return

        # Regular drag → Y-only zoom/pan
        # Skip tiny drags to prevent “zoom out on click”
        if ev.isFinish() and not self.dragging:
            ev.ignore()
            return

        self.dragging = True
        super().mouseDragEvent(ev, axis)

        if ev.isFinish():
            self.dragging = False
        else:
            self.dragging_signal.emit(True)
            

    def clear_selection(self):
        if self.lr is not None:
            self.removeItem(self.lr)
            self.lr = None
        self.shift_dragging = False

    def wheelEvent(self, ev, axis=None):
        # Y-only zoom unless Shift is pressed
        if ev.modifiers() == QtCore.Qt.ShiftModifier:
            ev.ignore()
            return
        factor = 0.9 if ev.delta() > 0 else 1.1
        y_min, y_max = self.viewRange()[1]
        center = 0.5 * (y_min + y_max)
        half = 0.5 * (y_max - y_min) * factor
        self.setYRange(center - half, center + half, padding=0)
        ev.accept()
        