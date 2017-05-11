import sys, os

import pymzml
import numpy as np
import pyqtgraph as pg
from PyQt4 import QtCore, QtGui
from pyqtgraph.Point import Point
from multiprocessing import Process, Queue

import mzML_viewer as MZML_V
import EIC_dialog as EIC_D
from xenophile.common import *


class EIC_dialog(QtGui.QDialog, EIC_D.Ui_Dialog):

    def __init__(self, files, parent = None):
        super(EIC_dialog, self).__init__(parent)
        #self.setStyleSheet('font-size: 10pt; font-family: Sans Serif;')
        self.setupUi(self)
        self.setWindowTitle('Plot EIC')

        self.make_dialog_connections()
        self.files = files
        self.msLevel_MS2.setChecked(False)

        # populate list widget with file names
        for i in range(len(files)):
            self.EIC_dialog_select_files.addItem(QtGui.QListWidgetItem(files[i].filename))
        return

    def make_dialog_connections(self):
        self.EIC_dialog_OK.clicked.connect(self.accept)
        self.EIC_dialog_cancel.clicked.connect(self.reject)
        return

    def returnValues(self):
        ''' get input from fields and return to mzML_view '''

        width = self.EIC_dialog_width.text()
        targets = [x.strip(',').strip() for x in str(self.EIC_dialog_targets.text()).split(',')]
        target_files = [str(x.text()) for x in self.EIC_dialog_select_files.selectedItems()]

        MS1 = self.msLevel_MS1.isChecked()
        MS2 = self.msLevel_MS2.isChecked()

        levels = []
        if MS1: levels.append(1)
        if MS2: levels.append(2)

        print levels
        return target_files, width, targets, levels

    # static method to create the dialog and return target_files, width, targets
    @staticmethod
    def get_EIC_specs(files, parent = None):
        dialog = EIC_dialog(files, parent)
        result = dialog.exec_()
        target_files, width, targets, levels = dialog.returnValues()
        return (target_files, width, targets, levels, result == QtGui.QDialog.Accepted)


class mzML_view (QtGui.QDialog, MZML_V.Ui_Dialog):

    def __init__(self, darkMode = True, parent = None):

        # Set PG plot background to white before GUI class init
        if not darkMode:
            pg.setConfigOption('background','w')
            #pg.setConfigOption('foreground', 'k')

        super(mzML_view, self).__init__(parent)
        #self.setStyleSheet('font-size: 10pt; font-family: Sans Serif;')
        self.setupUi(self)
        self.setWindowTitle('mzML File Explorer')

        ''' treewidget parameters '''
        self.list.setHeaderHidden(True)
        self.column = 0
        self.file_counter = 0
        self.parent = self.list

        ''' define plot details '''
        # add two plots to traces graphicsview object - EIC and MS
        self.EIC = self.traces.addPlot(row=1, col=0)
        self.MS = self.traces.addPlot(row=2, col=0)
        self.iline = pg.InfiniteLine(angle = 90, movable = True, pen = 'b')

        # EIC
        self.EIC.showGrid(x = True, y = True)
        self.EIC.showLabel('bottom', show = True)
        self.EIC.showLabel('left', show = True)
        self.EIC.showLabel('right', show = True)
        self.EIC.setLabel(axis = 'bottom', text = 'Retention Time (min)')
        self.EIC.setLabel(axis = 'left', text = 'Intensity')
        self.EIC.setLabel(axis = 'top', text = '')
        self.EIC.setLabel(axis = 'right', text = '')

        # MS
        self.MS.showGrid(x = True, y = True)
        self.MS.showLabel('bottom', show = True)
        self.MS.showLabel('left', show = True)
        self.MS.showLabel('right', show = True)
        self.MS.setLabel(axis = 'bottom', text = 'm/z')
        self.MS.setLabel(axis = 'left', text = 'Intensity')
        self.MS.setLabel(axis = 'top', text = '')
        self.MS.setLabel(axis = 'right', text = '')

        self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        self.nColors = len(self.colors)
        self.countColors = 0


        ''' class level variables '''
        self.files = [] # list of all file objects loaded
        self.fileCounter = 1
        self.iline_region_min = 0
        self.iline_region_max = 0
        self.iline_set = False
        self.iline_position = None
        self.active_trace = []
        self.stored_spectra = []

        ''' make connections '''
        self.make_mzML_connections()
        self.fileOpenDialogPath = os.path.expanduser('~')

        '''
        NOTES

        SO question on checking check status of treewidget items
        http://stackoverflow.com/questions/26963786/pyqt-get-list-of-all-checked-in-qtreewidget

        SO question on creating custom data entry dialog boxes
        http://stackoverflow.com/questions/5760622/pyqt4-create-a-custom-dialog-that-returns-parameters

        retrieving data from custom dialog
        http://stackoverflow.com/questions/18196799/how-can-i-show-a-pyqt-modal-dialog-and-get-data-out-of-its-controls-once-its-clo
        '''

    def make_mzML_connections(self):
        self.load_button.clicked.connect(self.load_data)
        self.list.itemChanged.connect(self.handle_changed)
        self.plot_chromatogram.clicked.connect(self.get_EICs)
        self.iline.sigPositionChangeFinished.connect(self.update_MS)
        self.iline.sigDragged.connect(self.update_MS)
        self.list.currentItemChanged.connect(self.get_active_trace)
        self.go_to_rt_button.clicked.connect(self.go_to_rt)
        self.closeHighlighted.clicked.connect(self.removeHighlightedFile)
        return

    def removeHighlightedFile(self):
        # get text of currently selected row in QTreeWidget

        # TODO -- need to destroy the relevant datafile objects
        # when the treeview item is removed
        root = self.list.invisibleRootItem()
        selected = self.list.selectedItems()
        if selected:
            listEntry = str(selected[0].text(0))
            self.files = [f for f in self.files if f.filename != listEntry]
        else:
            return
        for item in self.list.selectedItems():
            (item.parent() or root).removeChild(item)
        self.get_active_trace()
        self.update_plots()
        self.update_MS()
        return

    def load_data(self):

        filenames = [str(QtGui.QFileDialog.getOpenFileName(self,
                                                            'Select mzML File',
                                                            self.fileOpenDialogPath
                                                            ))]
        self.fileOpenDialogPath = os.path.dirname(filenames[0])
        for f in filenames:
            f = mzML_file(f, self.fileCounter)
            # spectra generator object for random access
            # NB - this only seems to work if no mzML file compressioin is performed
            f.spec = pymzml.run.Reader(f.filename)

            self.files.append(f)
            self.fileCounter += 1

        # extract initial TICs for selected files
        response = self.get_TIC()

        # add files to treewidget
        self.add_item(self.list)

        # update_plots
        self.update_plots()

        # check if this is the first chromatogram plotted
        if response is not None:
            min_rt, max_rt = response

            # calculate initial iline position
            position = min_rt + 0.01 * max_rt

            # add infinite line to EIC plot
            self.iline.setValue(position)
            self.iline_position = position
            self.EIC.addItem(self.iline, ignoreBounds=True)
        #self.list.setCurrentItem(0)
        return

    def get_TIC(self):
        '''
        plot EIC for each file in self.files list
            - adds self.RTs, and self.TIC np arrays to each mzML_file instance
            - self.RTs used for RT indexing of iLine position
        '''
        # create a list of lower bound RT values for setting linear region boundaries
        RT_bounds = []
        print 'Getting TIC'
        for filei in self.files:
            if not filei.new:
                continue
            # add TIC attributes
            filei.TIC_color = self.colors[self.countColors % self.nColors]
            self.countColors += 1

            msrun = pymzml.run.Reader(filei.filename)

            try:
                # try to get TIC data from mzML chromatogram - only seems to work if no mzML compression has been used
                tic = msrun['TIC']
                filei.RTs = np.asarray(tic.mz) # retention times are stored as .mz attribute for some reason
                filei.TIC = np.asarray(tic.i)

                # add boundary values to list
                RT_bounds.append(max(filei.RTs))
                RT_bounds.append(min(filei.RTs))

                continue # successful - move on to next file

            except Exception, e:
                print str(e)
                pass # build TIC from scratch

            # if here, mzML file does not have 'TIC' entry
            # assemble from scratch
            filei.RTs = np.empty(0)
            filei.TIC = np.empty(0)
            filei.levels = np.empty(0).astype(int) # array of indices of MSN levels

            del msrun
            spectra = readSpectra(filei.filename)

            for spectrum in spectra:


                time, mzs, ints, lvl = spectrum
                filei.RTs = np.append(filei.RTs, time)
                filei.TIC = np.append(filei.TIC, np.sum(ints))
                filei.levels = np.append(filei.levels, lvl)

            RT_bounds.append(np.amax(filei.RTs))
            RT_bounds.append(np.amin(filei.RTs))

            filei.ms1RTs = filei.RTs[np.where(filei.levels == 1)]
            filei.ms2RTs = filei.RTs[np.where(filei.levels == 2)]

            assert filei.levels.shape == filei.RTs.shape

        min_rt = min(RT_bounds)
        max_rt = max(RT_bounds)

        if self.iline_set == False:
            self.iline_set = True
            return [min_rt, max_rt]
        else:
            return None

    def get_active_trace(self):
        # get text of currently selected row in QTreeWidget
        item = self.list.currentItem()

        # get parent of this item
        try:
            item_text = str(item.text(0))
            parent = str(item.parent().text(0))
        except AttributeError:
            self.active_trace = []
            return

        # find matching file
        selectedFile = None
        for filei in self.files:
            if filei.filename == parent:
                selectedFile = filei
                break

        # update active trace
        del self.active_trace[:]

        if 'TIC' in item_text:
            trace_type = 'TIC'
            self.active_trace = [
                                trace_type,
                                selectedFile,
                                None
                                ]

        elif 'EIC' in item_text:
            for key, trace in filei.plot_data.iteritems():
                if key == item_text:
                    trace_type = 'EIC'
                    self.active_trace = [
                                        trace_type,
                                        selectedFile,
                                        trace
                                        ]

        return

    def get_stored_spectrum(self, n):

        for s in self.stored_spectra:
            if s.index == n:
                return s
        else:
            return None

    def update_MS(self):
        '''
        Called when iLine position is being moved
        - might bind left and right arrow keys to incrementally move iLine position
        '''
        if len(self.active_trace) == 0:
            return

        # unpack active_trace list
        Type, filei, trace = self.active_trace

        # clear MS plot
        self.MS.clear()

        # find boudaries of selected region
        iline_val = float(self.iline.value())

        if Type == 'TIC' or Type == 'EIC' and trace.levels == [1,2]:
            # find index of minimum absolute difference from iLine position and RT array
            index = np.argmin(np.absolute(filei.RTs - iline_val))

        elif Type == 'EIC':
            mask = (filei.levels == trace.levels[0])
            subset_id = np.argmin(np.absolute(filei.RTs[mask] - iline_val))
            index = np.arange(filei.RTs.shape[0])[mask][subset_id]

        # pull spectrum with this index
        try:
            spec = filei.spec[index+1]
            mz, i = self.zero_fill(np.asarray(spec.mz), np.asarray(spec.i))

        except: # mzML file not indexed

            # try to get spectrum from memory
            spectrum = self.get_stored_spectrum(index)

            if spectrum:
                mz, i = self.zero_fill(spectrum.mz, spectrum.i)
            else:
                # iterate through
                save_spectra = 100

                self.stored_spectra = [] # clear data
                spectra = self.readSpectra(filei.filename, None, stopIter = index + save_spectra)

                for n, spectrum in enumerate(spectra):
                    if n > index - save_spectra and n < index + save_spectra:
                        spectrum.index = n
                        self.stored_spectra.append(spectrum)
                    if n != index:
                        mz, i = self.zero_fill(spectrum.mz, spectrum.i)

        # plot spectrum
        self.MS.plot(x = mz, y = i, pen = 'b')

        pg.QtGui.QApplication.processEvents() # force complete plot redraw
        # or: app..processEvents() if running from top-level script
        return

    def zero_fill(self, xData, yData):
        x = np.repeat(xData, 3)
        y = np.dstack((np.zeros(yData.shape[0]), yData, np.zeros(yData.shape[0]))).flatten()
        return x, y

    def get_EICs(self):
        # call EIC_dialog class
        dialog = EIC_dialog(self.files)
        selected_files, width, targets, levels, result = dialog.get_EIC_specs(self.files)

        # result = boolean - True if return values exist, false if not
        if result:

            # build list of file objects from filenames in target_files list
            target_files = []
            for filei in self.files:
                if filei.filename in selected_files:
                    target_files.append(filei)

            # get EIC data
            self.plot_EIC(target_files, width, targets, levels)

            # add EIC element to treewidget child list
            for filei in target_files:
                for key, value in filei.plot_data.iteritems():

                    if value.isInTreeView == False:
                        self.addChild(filei.tree_parent, self.column, key, 'data Type A')
                        value.isInTreeView = True

            # update plots
            self.update_plots()

        else:
            return

    def go_to_rt(self):

        # move iLine ot provided rt
        target_rt = str(self.go_to_rt_field.text())
        self.iline.setValue(target_rt)

        self.update_MS()
        return

    def plot_EIC(self, target_files, EIC_width, targets, levels):

        print('Extracting EIC...')

        for filei in target_files:
            print('getting EICs for: %s' %filei.filename)

            EICs = []
            getLevels = False

            # check if mslevel array exists - this will False if the TIC was loaded from mzML chromatogram
            if not hasattr(filei, 'levels'):
                filei.levels = np.empty(filei.RTs.shape[0]+1).astype(int) # array of indices of MSN levels
                getLevels = True

            # prepare targeting data
            for target in targets:
                EIC = EIC_trace()
                EIC.type = 'EIC'
                EIC.levels = levels
                EIC.target = float(target)
                EIC.ref_string = 'EIC_%s' %target
                EIC.ll = float(target) - float(EIC_width)
                EIC.hl = float(target) + float(EIC_width)
                EIC.isInTreeView = False
                EIC.data = np.empty(0)
                EIC.color = self.colors[self.countColors % self.nColors]
                self.countColors += 1
                EICs.append(EIC)

            # instantiate spectrum generator
            #msrun = pymzml.run.Reader(filei.filename)
            count = 0

            spectra = self.readSpectra(filei.filename, None)
            # get data
            for n, spectrum in enumerate(spectra):

                time, mzs, ints, lvl = spectrum.rt, spectrum.mz, spectrum.i, spectrum.lvl

                if getLevels: filei.levels[n] = lvl

                # iterate through MS1 spectra
                if lvl in levels:

                    # add 0 entry to each EIC data array
                    for EIC in EICs:
                        points = np.where((mzs > EIC.ll) & (mzs < EIC.hl))
                        EIC.data = np.append(EIC.data, np.sum(ints[points]))

                count += 1

            if not hasattr(filei, 'ms1RTs'):
                filei.ms1RTs = filei.RTs[np.where(filei.levels == 1)]
                filei.ms2RTs = filei.RTs[np.where(filei.levels == 2)]

            print('number of spectra in MSrun is: %s' %count)
            # add results to file object
            for EIC in EICs:
                filei.plot_data[EIC.ref_string] = EIC

            print('Done!')
        return

    def add_item(self, parent):
        ''' for each file add a treewidget parent node '''
        for filei in self.files:
            if not filei.new: continue
            filei.new = False
            #filei.tree_parent = self.list.addTopLevelItem(QtGui.QTreeWidgetItem(parent, [filei.filename]))
            filei.tree_parent = self.addParent(self.parent, self.column, filei.filename, 'data Clients')
            self.addChild(filei.tree_parent, self.column, 'TIC', 'data Type A')
        return

    def addParent(self, parent, column, title, data):
        item = QtGui.QTreeWidgetItem(parent, [title])
        item.setData(column, QtCore.Qt.UserRole, data)
        item.setChildIndicatorPolicy(QtGui.QTreeWidgetItem.ShowIndicator)
        item.setExpanded (True)
        return item

    def addChild(self, parent, column, title, data):
        item = QtGui.QTreeWidgetItem(parent, [title])
        item.setData(column, QtCore.Qt.UserRole, data)
        item.setCheckState (column, QtCore.Qt.Checked)
        return item

    def handle_changed(self, item, column):
        #if item.checkState(column) == QtCore.Qt.Checked:
        #    print "checked", item, item.text(column)
        #f item.checkState(column) == QtCore.Qt.Unchecked:
        #    print "unchecked", item, item.text(column)

        self.update_plots()
        return

    def get_checked_items(self):

        # get list of all checked elements
        root = self.list.invisibleRootItem()
        file_count = root.childCount()
        checked_elements = []

        for i in range(file_count):
            filei = root.child(i)
            num_traces = filei.childCount()
            for n in range(num_traces):
                child = filei.child(n)
                if child.checkState(0) == QtCore.Qt.Checked:
                    checked_elements.append([str(filei.text(0)), str(child.text(0))])

        return checked_elements

    def update_plots(self):

        checked_elements = self.get_checked_items()

        # clear plots
        self.EIC.clear()

        # clearing EIC also erases iline -> replace this
        self.EIC.addItem(self.iline, ignoreBounds=True)

        try:
            self.eicLGD.scene().removeItem(self.eicLGD)
        except:
            pass

        self.eicLGD = self.EIC.addLegend()

        # draw new plots
        for e in checked_elements:
            for filei in self.files:
                if filei.filename == e[0]:
                    if e[1] != 'TIC':
                        for key, value in filei.plot_data.iteritems():
                            if key == e[1] and value.type == 'EIC':

                                if value.levels == [1]: rts = filei.ms1RTs
                                if value.levels == [2]: rts = filei.ms2RTs
                                if value.levels == [1,2]: rts = filei.RTs
                                self.EIC.plot(x = rts, y = value.data, pen = value.color, name = value.ref_string)
                    else:
                        self.EIC.plot(x = filei.RTs, y = filei.TIC, pen = filei.TIC_color, name = str(filei.filename.split('/')[-1]))

        return

    def readSpectra(self, mzml_file, msLevel, stopIter = None):

        import pymzml

        msrun = pymzml.run.Reader(str(mzml_file))
        for n, spectrum in enumerate(msrun):

            # only consider MS1 level
            if msLevel:
                if spectrum['ms level'] != msLevel: continue

            lvl = spectrum['ms level']

            try:
                time = spectrum['scan time']
            except:
                try:
                    time = spectrum['scan start time']
                except:
                    print 'skipping spectrum'
                    #print spectrum.keys()
                    continue

            try:

                if stopIter and n >= stopIter: break
                mzs = np.array(spectrum.mz, dtype="float32")
                ints = np.array(spectrum.i, dtype ='uint64')
                yield Spectrum(time, mzs, ints, lvl)

            except Exception, e:
                continue

class EIC_trace(object):
    def __init__(self):
        return

class mzML_file(object):
    def __init__(self, filename, i):
        self.filename = filename
        self.index = i
        self.plot_data = {}
        self.new = True
        return

class mzML_trace_data(object):
    def __init__(self):
        return

class Spectrum(object):
    def __init__(self, rt, mzs, ints, lvl):
        self.rt = rt
        self.mz = mzs
        self.i = ints
        self.lvl = lvl

def get_files():
        import qdarkstyle as qds

        dialog = QtGui.QFileDialog()
        dialog.setFileMode(QtGui.QFileDialog.AnyFile)
        dialog.setStyleSheet(qds.load_stylesheet(pyside = False))
        dialog.setViewMode(QtGui.QFileDialog.Detail)
        if dialog.exec_():
            filenames = dialog.selectedFiles()
        return filenames
