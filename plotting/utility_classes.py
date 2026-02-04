from PySide6 import QtWidgets
from PySide6 import QtCore
from PySide6.QtGui import QKeySequence, QShortcut
import pyqtgraph as pg
from datetime import datetime, timedelta
from obspy import UTCDateTime


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

    def clear_permanent_selection(self):
        if self.lr_permanent is not None:
            self.removeItem(self.lr_permanent)
            self.lr_permanent = None

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
    
class CatalogDialog(QtWidgets.QDialog):
    LOAD, SAVE, NEW, CANCEL = range(4)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Catalog")

        self.choice = self.CANCEL

        layout = QtWidgets.QVBoxLayout(self)

        label = QtWidgets.QLabel("What would you like to do?")
        layout.addWidget(label)

        btn_layout = QtWidgets.QHBoxLayout()

        load_btn = QtWidgets.QPushButton("Load Catalog")
        save_btn = QtWidgets.QPushButton("Save Catalog")
        new_btn  = QtWidgets.QPushButton("New Catalog")
        cancel_btn = QtWidgets.QPushButton("Cancel")

        load_btn.clicked.connect(lambda: self._done(self.LOAD))
        save_btn.clicked.connect(lambda: self._done(self.SAVE))
        new_btn.clicked.connect(lambda: self._done(self.NEW))
        cancel_btn.clicked.connect(lambda: self._done(self.CANCEL))

        for b in (load_btn, save_btn, new_btn, cancel_btn):
            btn_layout.addWidget(b)

        layout.addLayout(btn_layout)

    def _done(self, choice):
        self.choice = choice
        self.accept()


class EventNoteDialog(QtWidgets.QDialog):
    def __init__(self, start_time: UTCDateTime, end_time: UTCDateTime, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Add Event")
        self.resize(450, 280)

        start_time = start_time.datetime
        end_time = end_time.datetime

        layout = QtWidgets.QVBoxLayout(self)

        # Time group
        time_group = QtWidgets.QGroupBox("Event Time")
        time_layout = QtWidgets.QFormLayout(time_group)

        self.start_edit = QtWidgets.QDateTimeEdit(start_time)
        self.end_edit   = QtWidgets.QDateTimeEdit(end_time)

        for edit in (self.start_edit, self.end_edit):
            edit.setCalendarPopup(True)
            edit.setDisplayFormat("yyyy-MM-dd HH:mm:ss.zzz")
            edit.setTimeSpec(QtCore.Qt.LocalTime)
            edit.setMinimumWidth(200)

        time_layout.addRow("Start:", self.start_edit)
        time_layout.addRow("End:", self.end_edit)

        self.duration_label = QtWidgets.QLabel()
        time_layout.addRow("Duration:", self.duration_label)

        layout.addWidget(time_group)

        # Note entry
        layout.addWidget(QtWidgets.QLabel("Event note:"))

        self.text_edit = QtWidgets.QTextEdit()
        self.text_edit.setPlaceholderText("Enter a note for this event…")
        layout.addWidget(self.text_edit)

        # Buttons
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        )
        layout.addWidget(buttons)

        self.ok_button = buttons.button(QtWidgets.QDialogButtonBox.Ok)

        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        # Signals
        self.start_edit.dateTimeChanged.connect(self._validate)
        self.end_edit.dateTimeChanged.connect(self._validate)

        self._validate()
        self.text_edit.setFocus()

    def _validate(self):
        start = self.start_edit.dateTime()
        end = self.end_edit.dateTime()

        duration_ms = start.msecsTo(end)
        valid = duration_ms > 0

        self.ok_button.setEnabled(valid)

        if valid:
            seconds = duration_ms / 1000.0
            self.duration_label.setText(f"{seconds:.3f} s")
            self.duration_label.setStyleSheet("")
        else:
            self.duration_label.setText("Invalid (end ≤ start)")
            self.duration_label.setStyleSheet("color: red;")

    def qdatetime_to_utcdatetime(self, qdt):
        # Convert to ObsPy UTCDateTime
        return UTCDateTime(qdt.toString("yyyy-MM-ddTHH:mm:ss.zzz"))
    
    def start_time(self) -> UTCDateTime:
        print(self.start_edit.dateTime())
        return self.qdatetime_to_utcdatetime(self.start_edit.dateTime())

    def end_time(self) -> UTCDateTime:
        return self.qdatetime_to_utcdatetime(self.end_edit.dateTime())

    def note(self) -> str:
        return self.text_edit.toPlainText().strip()


class EventCatalogModel(QtCore.QAbstractTableModel):
    def __init__(self, catalog: dict, parent=None):
        super().__init__(parent)
        self._catalog = catalog
        self._keys = list(catalog.keys())  # row identifiers
        # Determine column names from first entry
        first_entry = next(iter(catalog.values()))
        self._columns = list(first_entry.keys())

    def rowCount(self, parent=None):
        return len(self._keys)

    def columnCount(self, parent=None):
        return len(self._columns)

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if not index.isValid():
            return None

        row = self._keys[index.row()]
        col = self._columns[index.column()]
        value = self._catalog[row].get(col, "")

        if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
            return str(value)
        return None

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role != QtCore.Qt.DisplayRole:
            return None
        if orientation == QtCore.Qt.Horizontal:
            return self._columns[section]
        else:
            # return string representation of the integer key
            return str(self._keys[section])

    def flags(self, index):
        if not index.isValid():
            return QtCore.Qt.ItemIsEnabled
        return QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        if role != QtCore.Qt.EditRole or not index.isValid():
            return False
        row = self._keys[index.row()]
        col = self._columns[index.column()]
        self._catalog[row][col] = value
        self.dataChanged.emit(index, index)
        return True

    # Optional: add/remove rows
    def insertRow(self, row, new_entry=None, parent=None):
        self.beginInsertRows(QtCore.QModelIndex(), row, row)
        key = len(self._keys)
        if new_entry is None:
            new_entry = {col: "" for col in self._columns}
        # set event_id if present
        if 'event_id' in self._columns:
            new_entry['event_id'] = key
        self._catalog[key] = new_entry
        self._keys.insert(row, key)
        self.endInsertRows()
        return True

    def removeRow(self, row, parent=None):
        self.beginRemoveRows(QtCore.QModelIndex(), row, row)
        key = self._keys.pop(row)
        del self._catalog[key]
        self.endRemoveRows()
        return True

    def moveRow(self, sourceRow, destinationRow):
        """Move a row from sourceRow to destinationRow"""
        if sourceRow == destinationRow or not (0 <= sourceRow < self.rowCount()) or not (0 <= destinationRow < self.rowCount()):
            return False

        self.beginMoveRows(QtCore.QModelIndex(), sourceRow, sourceRow, QtCore.QModelIndex(), destinationRow + (1 if destinationRow > sourceRow else 0))
        key = self._keys.pop(sourceRow)
        self._keys.insert(destinationRow, key)
        self.endMoveRows()
        return True
    
    def resetIndices(self, proxy_model: QtCore.QSortFilterProxyModel = None):
        """
        Renumber the integer keys sequentially according to current
        visible order (proxy model) or source model order if no proxy.
        """
        self.beginResetModel()

        if proxy_model is not None:
            # Get source row order according to proxy
            source_rows_in_order = [proxy_model.mapToSource(proxy_model.index(row, 0)).row()
                                    for row in range(proxy_model.rowCount())]
        else:
            # No proxy; default order
            source_rows_in_order = list(range(len(self._keys)))

        # Build new catalog with sequential keys
        new_keys = list(range(len(source_rows_in_order)))
        new_catalog = {}
        new_keys_map = {}

        for new_key, source_row in zip(new_keys, source_rows_in_order):
            old_key = self._keys[source_row]
            new_catalog[new_key] = self._catalog[old_key]
            new_keys_map[source_row] = new_key

        # Assign new catalog and keys
        self._catalog = new_catalog
        self._keys = new_keys

        self.endResetModel()

class NumericEventProxy(QtCore.QSortFilterProxyModel):
    def lessThan(self, left, right):
        l_data = left.data()
        r_data = right.data()
        print(f"Comparing '{l_data}' with '{r_data}'")
        try:
            # try to extract floating point numbers from strings
            l_num = float(l_data)
            r_num = float(r_data)
            return l_num < r_num
        except Exception:
            return str(l_data) < str(r_data)

class CatalogEditor(QtWidgets.QDialog):
    def __init__(self, catalog: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Event Catalog Editor")
        self.resize(900, 400)

        layout = QtWidgets.QVBoxLayout(self)

        self.model = EventCatalogModel(catalog)

        # Proxy model for sorting/filtering
        self.proxy_model = NumericEventProxy()
        self.proxy_model.setSourceModel(self.model)
        self.proxy_model.setSortCaseSensitivity(QtCore.Qt.CaseInsensitive)
        self.proxy_model.setDynamicSortFilter(True)

        # Table view uses the proxy model
        self.view = QtWidgets.QTableView()
        self.view.setModel(self.proxy_model)
        self.view.setSortingEnabled(True)
        self.proxy_model.sort(-1)
        self.view.resizeColumnsToContents()
        self.view.resizeRowsToContents()
        layout.addWidget(self.view)

        # Buttons
        buttons = QtWidgets.QHBoxLayout()
        add_btn = QtWidgets.QPushButton("Add Event")
        #remove_btn = QtWidgets.QPushButton("Remove Selected")
        up_btn = QtWidgets.QPushButton("Move Up")
        down_btn = QtWidgets.QPushButton("Move Down")
        reset_idx_btn = QtWidgets.QPushButton("Reset Indices")
        clear_sort_btn = QtWidgets.QPushButton("Clear Sorting")
        save_btn = QtWidgets.QPushButton("Save & Close")
        reject_btn = QtWidgets.QPushButton("Cancel")
        buttons.addWidget(add_btn)
        #buttons.addWidget(remove_btn)
        buttons.addWidget(clear_sort_btn)
        buttons.addWidget(up_btn)
        buttons.addWidget(down_btn)
        buttons.addWidget(reset_idx_btn)
        buttons.addStretch()
        buttons.addWidget(reject_btn)
        buttons.addWidget(save_btn)
        layout.addLayout(buttons)

        add_btn.clicked.connect(self.add_event)
        #remove_btn.clicked.connect(self.remove_event)
        clear_sort_btn.clicked.connect(self.clear_sorting)
        save_btn.clicked.connect(self.accept)
        up_btn.clicked.connect(lambda: self.move_selected(-1))
        down_btn.clicked.connect(lambda: self.move_selected(1))
        reset_idx_btn.clicked.connect(lambda: self.model.resetIndices(self.proxy_model))
        reject_btn.clicked.connect(self.reject)

        # bind delete key to remove_event
        delete_shortcut = QShortcut(QKeySequence.Delete, self.view)
        backspace_shortcut = QShortcut(QKeySequence.Backspace, self.view)
        delete_shortcut.activated.connect(self.remove_event)
        backspace_shortcut.activated.connect(self.remove_event)
    
    def clear_sorting(self):
        self.view.horizontalHeader().setSortIndicator(-1, QtCore.Qt.AscendingOrder)
        self.proxy_model.sort(-1)
    
    def is_sorted(self):
        return self.proxy_model.sortColumn() != -1

    def add_event(self):
        row = self.model.rowCount()
        self.model.insertRow(row)
        # Scroll to new row
        index = self.model.index(row, 0)
        self.view.scrollTo(index)
        self.view.edit(index)

    def remove_event(self):
        selection = self.view.selectionModel().selectedRows()
        if not selection:
            return

        # Map proxy rows to source rows
        source_rows = sorted([self.proxy_model.mapToSource(idx).row() for idx in selection], reverse=True)
        for row in source_rows:
            self.model.removeRow(row)

        
    def move_selected(self, direction):
        if self.is_sorted():
            QtWidgets.QMessageBox.information(
                self,
                "Reordering Disabled",
                "Clear sorting before manually reordering rows."
            )
            return

        selection = self.view.selectionModel().selectedRows()
        if not selection:
            return

        source_rows = [self.proxy_model.mapToSource(idx).row() for idx in selection]

        # Sort so moves don't collide
        source_rows.sort(reverse=(direction > 0))

        for row in source_rows:
            new_row = row + direction
            if 0 <= new_row < self.model.rowCount():
                self.model.moveRow(row, new_row)

        # Restore selection
        self.view.clearSelection()
        for row in source_rows:
            idx = self.proxy_model.mapFromSource(
                self.model.index(row + direction, 0)
            )
            self.view.selectRow(idx.row())


    def get_catalog(self) -> dict:
        # convert start_time and end_time back to UTCDateTime
        for event in self.model._catalog.values():
            if isinstance(event['start_time'], str):
                event['start_time'] = UTCDateTime(event['start_time'])
            if isinstance(event['end_time'], str):
                event['end_time'] = UTCDateTime(event['end_time'])
        return self.model._catalog

