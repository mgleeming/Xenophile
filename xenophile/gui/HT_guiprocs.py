import sys, os

import HT_search
import xenophile.libs.hitime.heatmap_subtraction as heatmap_subtraction
import xenophile.libs.hitime.HTS_resutls_file_parser as HTS_FP

from PyQt4 import QtCore, QtGui
import pyqtgraph as pg
from pyqtgraph.Point import Point

import numpy as np
from xenophile.common import *

class HT_search (QtGui.QDialog, HT_search.Ui_Dialog):

	def __init__(self, darkMode = True, parent = None, runner = None):

		# Set PG plot background to white before GUI class init
		if not darkMode:
			pg.setConfigOption('background','w')
			#pg.setConfigOption('foreground', 'k')

		super(HT_search, self).__init__(parent)

		self.setupUi(self)
		self.setWindowTitle('HiTIME Search')

		if runner:
			self.runner, self.runner_thread, self.q = runner

		self.timer = QtCore.QTimer()
		self.timer.timeout.connect(self.check_response)


		''' class level variables'''
		self.HT_search_list = []
		self.hits = [] # stores hitime hits for RV tab
		self.raw_data = None # raw hitime output data for RV tab
		self.raw_HT_RP_data = None # HT_RP raw data holding list
		self.refined_HT_RP_data = None # RP refined data

		self.RP_HT_file = None
		self.RP_mzML_file = None
		self.HT_RP_outputFileName = None

		self.ht_bs_treatment_data = None
		self.ht_bs_control_data = None
		self.ht_bs_output_data = None

		self.ht_bs_rtTolerance = 0.3
		self.ht_bs_mzTolerance = 0.2
		self.ht_bs_scoreCutoff = 0

		self.HT_BS_mz_tol.setText(str(self.ht_bs_mzTolerance))
		self.HT_BS_rt_tol.setText(str(self.ht_bs_rtTolerance))
		self.HT_BS_score_cutoff.setText(str(self.ht_bs_scoreCutoff))

		''' set plot config options '''

                self.makeRVPlots()

		''' set tablewidget geometries and selection behaviours '''
		self.HT_RV_hitlist.resizeRowsToContents()
		self.HT_RV_hitlist.resizeColumnsToContents()
		self.HT_RV_hitlist.horizontalHeader().setStretchLastSection(True)

		self.HTS_input_file_table.resizeRowsToContents()
		self.HTS_input_file_table.resizeColumnsToContents()
		self.HTS_input_file_table.horizontalHeader().setStretchLastSection(True)

		self.HTS_define_scan.resizeRowsToContents()
		self.HTS_define_scan.resizeColumnsToContents()
		self.HTS_define_scan.horizontalHeader().setStretchLastSection(True)

		self.HT_RV_accepted_list.resizeRowsToContents()
		self.HT_RV_accepted_list.resizeColumnsToContents()
		self.HT_RV_accepted_list.horizontalHeader().setStretchLastSection(True)

		''' set selection model for TableWidget items - want to select the entire row when a cell is clicked '''
		self.HT_RV_hitlist.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
		self.HT_RV_accepted_list.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)

		''' set tablewidget view - remove row counter column '''
		self.HT_RV_hitlist.verticalHeader().setVisible(False)
		self.HT_RV_accepted_list.verticalHeader().setVisible(False)

		''' set infinite lines for MS L/H markers '''
		self.loLine = pg.InfiniteLine(angle = 90, movable = False, pen = 'r')
		self.hiLine = pg.InfiniteLine(angle = 90, movable = False, pen = 'r')

		'''
		Plotting functios and data for HT postprocessing
		'''
		# HM
		self.RP_heat_map.showGrid(x = True, y = True)
		self.RP_heat_map.showLabel('bottom', show = True)
		self.RP_heat_map.showLabel('top', show = True)
		self.RP_heat_map.showLabel('left', show = True)
		self.RP_heat_map.showLabel('right', show = True)
		self.RP_heat_map.setLabel(axis = 'bottom', text = 'm/z')
		self.RP_heat_map.setLabel(axis = 'left', text = 'Retention Time (s)')
		self.RP_heat_map.setLabel(axis = 'top', text = '')
		self.RP_heat_map.setLabel(axis = 'right', text = '')

		# Histogram
		self.RP_histogram.showGrid(x = True, y = True)
		self.RP_histogram.showLabel('bottom', show = True)
		self.RP_histogram.showLabel('top', show = True)
		self.RP_histogram.showLabel('left', show = True)
		self.RP_histogram.showLabel('right', show = True)
		self.RP_histogram.setLabel(axis = 'bottom', text = 'Score')
		self.RP_histogram.setLabel(axis = 'left', text = 'Counts')
		self.RP_histogram.setLabel(axis = 'top', text = '')
		self.RP_histogram.setLabel(axis = 'right', text = '')

		'''
		Plot config options for subtraction tab
		'''
		# HM
		self.HT_BS_HM.showGrid(x = True, y = True)
		self.HT_BS_HM.showLabel('bottom', show = True)
		self.HT_BS_HM.showLabel('top', show = True)
		self.HT_BS_HM.showLabel('left', show = True)
		self.HT_BS_HM.showLabel('right', show = True)
		self.HT_BS_HM.setLabel(axis = 'bottom', text = 'm/z')
		self.HT_BS_HM.setLabel(axis = 'left', text = 'Retention Time (s)')
		self.HT_BS_HM.setLabel(axis = 'top', text = '')
		self.HT_BS_HM.setLabel(axis = 'right', text = '')

		''' form GUI connections '''
		self.HTS_gui_connections()
		self.HT_BS_connections()
		self.HT_RP_connections()
		self.HT_RV_connections()

		''' field activation defaults '''
		self.EICtoggled()
		self.fileOpenDialogPath = os.path.expanduser('~')

	def HTS_gui_connections(self):
		self.HTS_load_input_files_button.clicked.connect(self.HTS_select_file)
		self.HTS_run_HT_search.clicked.connect(self.HTS_run_HiTIME)
		self.HTS_apply_params.clicked.connect(self.HTS_update_search_params)
		self.HTS_input_file_table.itemSelectionChanged.connect(self.HTS_update_define_search)
		self.HTS_set_output_file.clicked.connect(self.HTS_set_output_directory)
		self.HTS_restore_defaults.clicked.connect(self.HTS_restore_default_params)
		self.HTS_remove_row_button.clicked.connect(self.HTS_remove_input_file)
		return

	def HTS_remove_input_file(self):
		highlighted_row = self.HTS_input_file_table.selectionModel().selectedRows()[0].row()
		selected_file = str(self.HTS_input_file_table.item(highlighted_row, 0).text())
		self.HT_search_list = [x for x in self.HT_search_list if x['inputFile'] != selected_file]
		self.HTS_input_file_table.removeRow(int(highlighted_row))
		try:
			self.HTS_input_file_table.selectRow(0)
		except:
			pass
		return

	def HTS_select_file(self):
		#self.lineEdit_2.setText(QtGui.QFileDialog.getOpenFileName())
		HT_file = [QtGui.QFileDialog.getOpenFileName(self, 'Select HiTIME File', self.fileOpenDialogPath)]
		self.fileOpenDialogPath = os.path.dirname(str(HT_file))

		defaultmzDelta = 6.0201
		defaultIR = 1
		defaultIC = 0
		defaultmzWidth = 150
		defaultrtWidth = 17
		defaultmzSigma = 1.5
		defaultrtSigma = 1.5
		defaultminSample = defaultrtWidth * defaultrtSigma/2.355

		# get row coung
		rc = self.HTS_input_file_table.rowCount()
		self.HTS_input_file_table.insertRow(rc)

		for i, f in enumerate(HT_file):
			f = str(f)

			file_entry = {
				'inputFile' : f,
				'threads' : 1,
				'outputFile' : os.path.splitext(f)[0] + '_out.dat',
				'mzDelta' : defaultmzDelta,
				'intensityRatio' : defaultIR,
				'mzWidth' : defaultmzWidth,
				'rtWidth' : defaultrtWidth,
				'rtSigma' : defaultrtSigma,
				'mzSigma' : defaultmzSigma,
				'minSample' : defaultminSample,
				'ppm' : 4,
				'logFile' : False,
				'noScore' : False,
				'removeLow' : False,
				'outDir' : None,
				'format' : 'mzml'
			}

			self.HT_search_list.append(file_entry)
			self.HTS_input_file_table.setItem(i + rc, 0, QtGui.QTableWidgetItem(f))

		self.HTS_input_file_table.selectRow(rc)
		self.HTS_input_file_table.resizeRowsToContents()
		self.HTS_input_file_table.resizeColumnsToContents()
		self.HTS_input_file_table.horizontalHeader().setStretchLastSection(True)

		self.HTS_update_define_search()
		return

	def HTS_clear_define_search(self):
		self.HTS_define_scan.setItem(0,0, QtGui.QTableWidgetItem(str('')))
		self.HTS_define_scan.setItem(0,1, QtGui.QTableWidgetItem(str('')))
		self.HTS_define_scan.setItem(0,2, QtGui.QTableWidgetItem(str('')))
		self.HTS_define_scan.setItem(0,3, QtGui.QTableWidgetItem(str('')))
		self.HTS_define_scan.setItem(0,4, QtGui.QTableWidgetItem(str('')))
		self.HTS_define_scan.setItem(0,5, QtGui.QTableWidgetItem(str('')))
		self.HTS_define_scan.setItem(0,6, QtGui.QTableWidgetItem(str('')))
		self.HTS_define_scan.setItem(0,7, QtGui.QTableWidgetItem(str('')))
		return

	def HTS_update_define_search(self):
		try:
			highlighted_row = self.HTS_input_file_table.selectionModel().selectedRows()[0].row()
			selected_file = str(self.HTS_input_file_table.item(highlighted_row, 0).text())
		except:
			self.HTS_clear_define_search()
			return

		search_params = None
		for entry in self.HT_search_list:
			if entry['inputFile'] == selected_file:
				search_params = entry
				break

		if not search_params: return
		#print search_params['inputFile']
		self.HTS_define_scan.setItem(0,0, QtGui.QTableWidgetItem(str(search_params['mzDelta'])))
		self.HTS_define_scan.setItem(0,1, QtGui.QTableWidgetItem(str(search_params['intensityRatio'])))
		self.HTS_define_scan.setItem(0,2, QtGui.QTableWidgetItem(str(search_params['mzWidth'])))
		self.HTS_define_scan.setItem(0,3, QtGui.QTableWidgetItem(str(search_params['rtWidth'])))
		self.HTS_define_scan.setItem(0,4, QtGui.QTableWidgetItem(str(search_params['mzSigma'])))
		self.HTS_define_scan.setItem(0,5, QtGui.QTableWidgetItem(str(search_params['rtSigma'])))
		self.HTS_define_scan.setItem(0,6, QtGui.QTableWidgetItem(str(search_params['threads'])))
		self.HTS_define_scan.setItem(0,7, QtGui.QTableWidgetItem(str(search_params['outputFile'])))

		self.HTS_define_scan.item(0,0).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,1).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,2).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,3).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,4).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,5).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,6).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,7).setTextAlignment(QtCore.Qt.AlignCenter)

		self.HTS_define_scan.resizeRowsToContents()
		self.HTS_define_scan.resizeColumnsToContents()
		self.HTS_define_scan.horizontalHeader().setStretchLastSection(True)
		return

	def HTS_update_search_params(self):
		highlighted_row = self.HTS_input_file_table.selectionModel().selectedRows()[0].row()
		selected_file = str(self.HTS_input_file_table.item(highlighted_row, 0).text())
		search_params = None
		for entry in self.HT_search_list:
			if entry['inputFile'] == selected_file:
				search_params = entry
				break
		if not search_params: return
		search_params['mzDelta'] = float(self.HTS_define_scan.item(0, 0).text())
		search_params['intensityRatio'] = float(self.HTS_define_scan.item(0, 1).text())
		search_params['mzWidth'] = float(self.HTS_define_scan.item(0, 2).text())
		search_params['rtWidth'] = float(self.HTS_define_scan.item(0, 3).text())
		search_params['mzSigma'] = float(self.HTS_define_scan.item(0, 4).text())
		search_params['rtSigma'] = float(self.HTS_define_scan.item(0, 5).text())
		search_params['threads'] = str(self.HTS_define_scan.item(0, 6).text())
		search_params['outputFile'] = str(self.HTS_define_scan.item(0, 7).text())
		return

	def HTS_set_output_directory(self):
		output_file = QtGui.QFileDialog.getSaveFileName(self, 'Select Output File', self.fileOpenDialogPath)
		self.fileOpenDialogPath = os.path.dirname(str(output_file))
		self.HTS_define_scan.setItem(0,7, QtGui.QTableWidgetItem(str(output_file)))
		return

	def HTS_restore_default_params(self):
		defaultmzDelta = 6.0201
		defaultIR = 1
		defaultIC = 0
		defaultmzWidth = 150
		defaultrtWidth = 17
		defaultmzSigma = 1.5
		defaultrtSigma = 1.5
		defaultminSample = defaultrtWidth * defaultrtSigma/2.355

		highlighted_row = self.HTS_input_file_table.selectionModel().selectedRows()[0].row()
		selected_file = str(self.HTS_input_file_table.item(highlighted_row, 0).text())
		outFile = os.path.splitext(selected_file)[0] + '_out.dat'

		self.HTS_define_scan.setItem(0,0, QtGui.QTableWidgetItem(str(defaultmzDelta)))
		self.HTS_define_scan.setItem(0,1, QtGui.QTableWidgetItem(str(defaultIR)))
		self.HTS_define_scan.setItem(0,2, QtGui.QTableWidgetItem(str(defaultmzWidth)))
		self.HTS_define_scan.setItem(0,3, QtGui.QTableWidgetItem(str(defaultrtWidth)))
		self.HTS_define_scan.setItem(0,4, QtGui.QTableWidgetItem(str(defaultmzSigma)))
		self.HTS_define_scan.setItem(0,5, QtGui.QTableWidgetItem(str(defaultrtSigma)))
		self.HTS_define_scan.setItem(0,6, QtGui.QTableWidgetItem(str(1)))
		self.HTS_define_scan.setItem(0,7, QtGui.QTableWidgetItem(str(outFile)))

		self.HTS_define_scan.item(0,0).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,1).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,2).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,3).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,4).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,5).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,6).setTextAlignment(QtCore.Qt.AlignCenter)
		self.HTS_define_scan.item(0,7).setTextAlignment(QtCore.Qt.AlignCenter)
		return

	def textEditAppend(self, f):
		self.HTS_status_text.appendPlainText(QtCore.QString(f))
		return

	def HTS_run_HiTIME(self):
		import xenophile.libs.hitime.hitime_methods as HTM
		run_job(self.runner, self.runner_thread, self.q, launch_HT_search, self.HT_search_list)
		self.timer.start(1000)
		return

	def check_response(self):

		while not self.q.empty():
			update = self.q.get()
			if update == 'done':
				time.sleep(0.1)
				self.runner.p.terminate()
				self.timer.stop()
				self.runner_thread.exit()
				print('\nResponse returned to guiprocs')
			else:
				self.textEditAppend(update)
				print update
		return

	'''
	HEATMAP SUBTRACTION TAB
	'''
	def HT_BS_connections(self):
		self.HT_BS_treatment_browse.clicked.connect(self.ht_bs_selectTreatment)
		self.HT_BS_control_browse.clicked.connect(self.ht_bs_selectControl)
		self.HT_BS_output_browse.clicked.connect(self.ht_bs_selectOutput)
		self.HT_BS_run.clicked.connect(self.ht_bs_runSubtraction)
		self.HT_BS_active_treatment.toggled.connect(self.ht_bs_updatePlot)
		self.HT_BS_active_control.toggled.connect(self.ht_bs_updatePlot)
		self.HT_BS_active_output.toggled.connect(self.ht_bs_updatePlot)
		return

	def ht_bs_selectTreatment(self):
		self.HT_BS_treatment_file = QtGui.QFileDialog.getOpenFileName(self, 'Select Treatment HiTIME File', self.fileOpenDialogPath)
		self.fileOpenDialogPath = os.path.dirname(str(self.HT_BS_treatment_file))
		self.HT_BS_treatment_field.setText(str(self.HT_BS_treatment_file))
		return

	def ht_bs_selectControl(self):
		self.HT_BS_control_file = QtGui.QFileDialog.getOpenFileName(self, 'Select Control HiTIME File', self.fileOpenDialogPath)
		self.fileOpenDialogPath = os.path.dirname(str(self.HT_BS_control_file))
		self.HT_BS_control_field.setText(str(self.HT_BS_control_file))
		return

	def ht_bs_selectOutput(self):
		self.HT_BS_outputFileName = str(QtGui.QFileDialog.getSaveFileName(self, 'Select Output File', self.fileOpenDialogPath))
		self.HT_BS_output_field.setText(self.HT_BS_outputFileName)
		return

	def ht_bs_updatePlot(self):

		filename = None

		# get checked radio button
		buttons = [self.HT_BS_active_treatment, self.HT_BS_active_control, self.HT_BS_active_output]
		filenames = [str(self.HT_BS_treatment_field.text()), str(self.HT_BS_control_field.text()), str(self.HT_BS_output_field.text())]
		for i, button in enumerate(buttons):
			if button.isChecked():
				# plot relevant data if available
				filename = filenames[i]
				break

		# exit if no file found
		if not filename: return

		self.HT_BS_HM.clear()

		htHeaders, data = read_hitime_files(
							filename,
							returnNp = True,
							retrieveTop = 10000,
							)
		self.HT_BS_HM.plot(x = data['mz'], y = data['rt'], pen = None, symbol = 'o')

		return

	def ht_bs_runSubtraction(self):
		# pull parameters

		rtTol = float(self.HT_BS_rt_tol.text())
		mzTol = float(self.HT_BS_mz_tol.text())
		scoreCutoff = float(self.HT_BS_score_cutoff.text())

		args = {
			'inTreatment' : self.HT_BS_treatment_file,
			'inControl' : self.HT_BS_control_file,
			'outFile' : self.HT_BS_outputFileName,
			'rtTolerance' : rtTol,
			'mzTolerance' : mzTol,
			'scoreCutoff' : scoreCutoff
		}

		run_job(self.runner, self.runner_thread, self.q, heatmap_subtraction.gui_init, args)
		self.timer.start(1000)
		return

	'''
	DATA POSTPROCESSING FUNCTIONS
	'''
	def HT_RP_connections(self):
		self.RP_HT_input_button.clicked.connect(self.parse_HT_data)
		self.RP_mascot_input_button.clicked.connect(self.select_mzml_file)
		self.RP_mzWidth.textChanged.connect(self.refine_data_plots)
		self.RP_rtWidth.textChanged.connect(self.refine_data_plots)
		self.HT_RP_rt_exclusion.textChanged.connect(self.refine_data_plots)
		self.RP_min_HT_score.textChanged.connect(self.refine_data_plots)
		self.RP_output_file_button.clicked.connect(self.set_output_file)
		self.RP_run_button.clicked.connect(self.run_data_postprocessing)
		self.HT_RP_plot_EICs.toggled.connect(self.EICtoggled)
		return

	def EICtoggled(self):
		eicState = self.HT_RP_plot_EICs.isChecked()
		self.RP_mzDelta.setEnabled(eicState)
		self.RP_EIC_width.setEnabled(eicState)
		self.RP_mascot_input_field.setEnabled(eicState)
		self.RP_mascot_input_button.setEnabled(eicState)
		return

	def set_output_file(self):
		self.HT_RP_outputFileName = str(QtGui.QFileDialog.getSaveFileName(self, 'Select Output File Name', self.fileOpenDialogPath))
		self.fileOpenDialogPath = os.path.dirname(str(self.HT_RP_outputFileName))
		self.RP_output_file_field.setText(self.HT_RP_outputFileName)
		return

	def select_mzml_file(self):
		self.RP_mzML_file = QtGui.QFileDialog.getOpenFileName(self, 'Select mzML Data File', self.fileOpenDialogPath)
		self.fileOpenDialogPath = os.path.dirname(str(self.RP_mzML_file))
		self.RP_mascot_input_field.setText(str(self.RP_mzML_file))
		return

	def parse_HT_data(self):

		'''
		Select HT input:
			Valid formats are:
				1) rt, mz, amp, score
				2) rt, mz, score

		Parse HT data into list of hitime_hit instances
			Attributes are
				.rt
				.mz
				.score
				.amp

		Bind data and plot histogram and heatmap

		Initially, all data is plotted which can then be refined by the user
		'''

		# get input file
		self.RP_HT_file = QtGui.QFileDialog.getOpenFileName(self, 'Select HiTIME Results File', self.fileOpenDialogPath)
		self.fileOpenDialogPath = os.path.dirname(str(self.RP_HT_file))
		self.RP_HT_input_field.setText(str(self.RP_HT_file))

		# run file parser
		self.HT_RP_headesr, self.raw_HT_RP_data = read_hitime_files(self.RP_HT_file)

		if not self.raw_HT_RP_data: return

		# parse data and plot histogram
		self.plot_score_histogram(self.raw_HT_RP_data)
		return

	def refine_data_plots(self):

		minScore = str(self.RP_min_HT_score.text())
		Drt = str(self.RP_rtWidth.text())
		Dmz = str(self.RP_mzWidth.text())
		rtExclusion = str(self.HT_RP_rt_exclusion.text())

		minScore = float(minScore) if minScore != '' else 0
		Drt = float(Drt) if Drt != '' else 0
		Dmz = float(Dmz) if Dmz != '' else 0
		rtExclusion = float(rtExclusion if rtExclusion != '' else 0.00)

		self.refined_HT_RP_data = get_HT_regions(self.raw_HT_RP_data, Drt, Dmz, minScore = minScore, rtExclusion = rtExclusion)

		self.plot_score_histogram(self.refined_HT_RP_data)

		return

	def plot_score_histogram(self, HT_data):

		# TODO >> this is a mess > fixme please

		# clear plots
		self.RP_histogram.clear()
		self.RP_heat_map.clear()

		# get a list of scores only
		scores = [float(x.score) for x in HT_data]

		max_score, min_score = max(scores), min(scores)

		# bin data - create 50 bins b/w min and max vals
		bin_counts, score_bins = np.histogram(scores, bins = np.linspace(min_score, max_score, 50))

		# plot data
		self.RP_histogram.plot(score_bins, bin_counts,stepMode = True, fillLevel = 0, brush = 'b')

		'''
		Hit Heatmap
		Note - increase plotting speed > only plot highest scoring 10,000 points
		'''
		rts = [float(x.rt) for x in HT_data[:10000]]
		mzs = [float(x.mz) for x in HT_data[:10000]]

		self.RP_heat_map.plot(x = mzs, y = rts, pen = None, symbol = 'o')

		# get histogram plot range limits
		histogram_x_range = (
					min_score,
					max_score
			)

		histogram_y_range = (
					0,
					np.log10(max(bin_counts))
			)

		# set plot ranges
		self.RP_histogram.setRange(xRange = histogram_x_range, yRange = histogram_y_range)
		# Set y axis of histogram to log scale
		self.RP_histogram.setLogMode(y = True)

		return

	def run_data_postprocessing(self):

		minScore = str(self.RP_min_HT_score.text())
		Drt = str(self.RP_rtWidth.text())
		Dmz = str(self.RP_mzWidth.text())
		mzDelta = str(self.RP_mzDelta.text())
		eicWidth = str(self.RP_EIC_width.text())
		rtExclusion = str(self.HT_RP_rt_exclusion.text())

		minScore = float(minScore) if minScore != '' else 0
		Drt = float(Drt) if Drt != '' else 0
		Dmz = float(Dmz) if Dmz != '' else 0
		mzDelta = float(mzDelta) if mzDelta != '' else 0
		eicWidth = float(eicWidth) if eicWidth != '' else 0.03
		rtExclusion = float(rtExclusion if rtExclusion != '' else 0.00)

		args = {
			'htIn' : str(self.RP_HT_file),
			'mzmlFile' : str(self.RP_mzML_file),
			'mzDelta' : mzDelta,
			'eicWidth' : eicWidth,
			'scoreCutoff' : minScore,
			'outFile' : str(self.HT_RP_outputFileName),
			'mzWidth' : Dmz,
			'rtWidth' : Drt,
			'peakList' : False,
			'rtExclusion' : rtExclusion,
			'plotEICs' : self.HT_RP_plot_EICs.isChecked(),
			'usePeptideIsotopeScaling' : self.HT_RP_use_peptide_isotope_scaling.isChecked()
		}

		import xenophile.libs.hitime.HT_search_postprocessing as HSP
		run_job(self.runner, self.runner_thread, self.q, HSP.guiRun, args)
		self.timer.start(1000)
		return

	'''
	CONNECTIONS AND FUNCTIONS FOR HITIME RESULTS TAB
	'''
	def HT_RV_connections(self):
		self.HT_RV_load_results.clicked.connect(self.parse_HT_results_file)
		self.HT_RV_hitlist.itemSelectionChanged.connect(lambda: self.HT_RV_update_plots(None))
		self.HT_RV_accepted_list.itemSelectionChanged.connect(lambda: self.HT_RV_update_plots(None, acceptedHit = True))
		self.HT_RV_up.clicked.connect(lambda: self.HT_RV_update_plots('up'))
		self.HT_RV_down.clicked.connect(lambda: self.HT_RV_update_plots('down'))
		self.HT_RV_reject.clicked.connect(lambda: self.HT_RV_update_plots('down'))
		self.HT_RV_accept.clicked.connect(self.add_hit_to_accepted_list)
		self.HT_RV_accept.clicked.connect(lambda: self.HT_RV_update_plots('down'))
		self.HT_RV_remove.clicked.connect(self.remove_hit_from_accepted_list)
		self.HT_RV_write.clicked.connect(self.write_accepted_to_file)
		self.RV_reset.clicked.connect(self.reset)
		self.HT_RV_HM.sigClicked.connect(self.pointClicked)
		return

	def pointClicked(self, plot, points):

		# find hit list index of clicked point
		for p in points:
			x = str(points[0]._data[0])
			y = str(points[0]._data[1])
			for i, h in enumerate(self.hits):
				if x == str(h.mz) and y == str(h.rt):
		 			# relevant result found - highlight row
		 			self.HT_RV_hitlist.selectRow(i)
					# update plots
					#self.HT_RV_update_plots(None)
					break
		return

	def keyPressEvent(self, event):
		# update plot ranges when the 'r' key is pressed

		keystroke = str(event.text())
		if keystroke == 'r' or keystroke == 'R':
			self.HT_RV_update_plots(None)
			return
		elif keystroke == '':
			self.close()
		else:
			return

	def parse_HT_results_file(self):

		inFile = QtGui.QFileDialog.getOpenFileName(self,
								'Specify Postprocessed HiTIME Reults File',
								self.fileOpenDialogPath
								)

		# safe browser dirname to defaults
		self.fileOpenDialogPath = os.path.dirname(str(inFile))

		# run file parser
		raw_HT_data, HT_hits, self.headers = HTS_FP.main(inFile)

		# store results in class level variables
		self.raw_data = raw_HT_data
		self.hits = [x for x in HT_hits]
		del raw_HT_data, HT_hits

		self.mzDelta = self.headers['mzDelta']

		# populate resutls table
		for i, hit in enumerate(self.hits):
			self.HT_RV_hitlist.insertRow(i)
			self.HT_RV_hitlist.setItem(i, 0, QtGui.QTableWidgetItem(str(hit.index)))
			self.HT_RV_hitlist.setItem(i, 1, QtGui.QTableWidgetItem(str(hit.rt)))
			self.HT_RV_hitlist.setItem(i, 2, QtGui.QTableWidgetItem('%.4f' %float(hit.mz)))
			self.HT_RV_hitlist.setItem(i, 3, QtGui.QTableWidgetItem('%.1f' %float(hit.score)))

			# centre text in tablewidget columns
			self.HT_RV_hitlist.item(i,0).setTextAlignment(QtCore.Qt.AlignCenter)
			self.HT_RV_hitlist.item(i,1).setTextAlignment(QtCore.Qt.AlignCenter)
			self.HT_RV_hitlist.item(i,2).setTextAlignment(QtCore.Qt.AlignCenter)
			self.HT_RV_hitlist.item(i,3).setTextAlignment(QtCore.Qt.AlignCenter)

		self.HT_HMX = [float(_.mz) for _ in self.hits]
		self.HT_HMY = [float(_.rt) for _ in self.hits]

		# set first row as highlighted
		self.HT_RV_hitlist.selectRow(0) #
		# NOTE: calling selectRow on the tableWidget item also calls the linked updatePlot method
		# ---> separately calling update_plots is unnecessary
		# # run update_plots to plot EIC and MS
		# self.HT_RV_update_plots(None)

		# update table measurements to fit data
		self.HT_RV_hitlist.resizeRowsToContents()
		self.HT_RV_hitlist.resizeColumnsToContents()
		self.HT_RV_hitlist.horizontalHeader().setStretchLastSection(True)
		return

	def HT_RV_update_plots(self, direction, acceptedHit = False):
		'''
		Update data analysis plots and stats when
			1) up or down button is clicked, or
			2) a new peptide entry is clicked in the results list
			3) MS spectrum type is changed using radio buttons
		'''

		if not acceptedHit:

			# remove any highlighting in accepted list
			self.HT_RV_accepted_list.clearSelection()

			# get selected row from RV_peptides table
			if direction is not None:
				index = self.HT_RV_hitlist.selectionModel().selectedRows()[0].row()
				if direction == 'up':
					self.HT_RV_hitlist.selectRow(index-1)
				elif direction == 'down':
					self.HT_RV_hitlist.selectRow(index+1)

			highlighted_row = self.HT_RV_hitlist.selectionModel().selectedRows()[0].row()
			hitNumber = int(str(self.HT_RV_hitlist.item(highlighted_row, 0).text()))

		else:
                        try:
                            # item selected in accepted list
			    highlighted_row = self.HT_RV_accepted_list.selectionModel().selectedRows()[0].row()
			    hitNumber = int(str(self.HT_RV_accepted_list.item(highlighted_row, 0).text()))
                        except:
                            return

		# get matching hit
		hit = None
		for x in self.hits:
			if x.index == hitNumber:
				hit = x
				break

		# clear existing plots
		self.HT_RV_MS.clear()
		self.HT_RV_EIC.clear()

		MS_mz = np.asarray(hit.MS_mz)
		MS_int = np.asarray(hit.MS_int)

		EIC_RT = np.asarray(hit.EIC_RT, dtype = 'float32')
		EIC_int_light = np.asarray(hit.EIC_int_light, dtype = 'uint64')
		EIC_int_heavy = np.asarray(hit.EIC_int_heavy, dtype = 'uint64')

		# set ranges
		xRangeMS = self.getXRange(hit.mz, 10, rangeType = 'fixed')
		xRangeRT = self.getXRange(hit.rt, 10, rangeType = 'pc')

		yRangeMS = self.getYRange(MS_mz, MS_int, xRangeMS, 'Processing SPECTRUM')

		yRangeRTL = self.getYRange(EIC_RT, EIC_int_light, xRangeRT, 'Processing LIGHT EIC')
		yRangeRTH = self.getYRange(EIC_RT, EIC_int_heavy, xRangeRT, 'Processing HEAVY EIC')

		yRangeRT = max([yRangeRTL, yRangeRTH])

		yRangeMS = ( 0 - yRangeMS * 0.05, yRangeMS)
		yRangeRT = ( 0 - yRangeRT * 0.05, yRangeRT)

		self.HT_RV_EIC.setLimits(xMin = float(min(EIC_RT)), xMax = float(max(EIC_RT)), yMin = yRangeMS[0], yMax = yRangeRT[1]*1.5)
		self.HT_RV_MS.setLimits(xMin = float(min(MS_mz)), xMax = float(max(MS_mz)), yMin = yRangeRT[0], yMax = yRangeMS[1]*1.5)

		self.HT_RV_MS.setRange(xRange = xRangeMS, yRange = yRangeMS)
		self.HT_RV_EIC.setRange(xRange = xRangeRT, yRange = yRangeRT)

		# set iLine positions and add to MS
		self.loLine.setValue(hit.mz)
		self.hiLine.setValue(hit.mz + self.mzDelta)

		self.HT_RV_MS.addItem(self.loLine, ignoreBounds = True)
		self.HT_RV_MS.addItem(self.hiLine, ignoreBounds = True)

		# plot MS
		self.HT_RV_MS.plot(
						x = np.repeat(MS_mz, 3),
						y = np.dstack((np.zeros(MS_int.shape[0]), MS_int, np.zeros(MS_int.shape[0]))).flatten(),
						pen = 'b',
						connect = 'all'
						)

		# plot EIC
		self.HT_RV_EIC.plot(x = EIC_RT, y = EIC_int_light, pen = 'b')
		self.HT_RV_EIC.plot(x = EIC_RT, y = EIC_int_heavy, pen = 'r')

		# plot heatmap
		self.HT_RV_HM.clear()

                self.HT_RV_HM.setData(x = self.HT_HMX, y = self.HT_HMY, pen = None, symbol = 'o', brush = 'b')
		self.HT_RV_HM.addPoints(x = [float(hit.mz)], y = [float(hit.rt)], symbol = 'o', brush = 'r')

		return

	def add_hit_to_accepted_list(self):
		''' take highlighted row from hit list and add to accepted hits '''

		# find highlighted row
		highlighted_row = self.HT_RV_hitlist.selectionModel().selectedRows()[0].row()
		hitNumber = int(self.HT_RV_hitlist.item(highlighted_row, 0).text())

		for hit in self.hits:
			if hit.index == hitNumber:
				rowCount = self.HT_RV_accepted_list.rowCount()
				self.HT_RV_accepted_list.insertRow(rowCount)
				self.HT_RV_accepted_list.setItem(rowCount, 0, QtGui.QTableWidgetItem(str(hit.index)))
				self.HT_RV_accepted_list.setItem(rowCount, 1, QtGui.QTableWidgetItem(str(hit.rt)))
				self.HT_RV_accepted_list.setItem(rowCount, 2, QtGui.QTableWidgetItem('%.4f' %float(hit.mz)))
				self.HT_RV_accepted_list.setItem(rowCount, 3, QtGui.QTableWidgetItem('%.1f' %float(hit.score)))

				# centre text in tablewidget columns
				self.HT_RV_accepted_list.item(rowCount, 0).setTextAlignment(QtCore.Qt.AlignCenter)
				self.HT_RV_accepted_list.item(rowCount, 1).setTextAlignment(QtCore.Qt.AlignCenter)
				self.HT_RV_accepted_list.item(rowCount, 2).setTextAlignment(QtCore.Qt.AlignCenter)
				self.HT_RV_accepted_list.item(rowCount, 3).setTextAlignment(QtCore.Qt.AlignCenter)
				break

		# update table measurements to fit data
		self.HT_RV_accepted_list.resizeRowsToContents()
		self.HT_RV_accepted_list.resizeColumnsToContents()
		self.HT_RV_accepted_list.horizontalHeader().setStretchLastSection(True)
		return

	def remove_hit_from_accepted_list(self):
		# find highlighted row
		highlighted_row = self.HT_RV_accepted_list.selectionModel().selectedRows()[0].row()
		self.HT_RV_accepted_list.removeRow(highlighted_row)
		return

	def write_accepted_to_file(self):
		# get output file name and location
		ofname = str(QtGui.QFileDialog.getSaveFileName(self, 'Set Output File', self.fileOpenDialogPath))

		# write data
		of1 = open(ofname, 'wt')

		of1.write('#- Validated HiTIME hit list:\n')

		# write postprocessing headers to results file
		of1.write('#- Postprocessing parameters\n')
		for k, v in self.headers.iteritems():
			of1.write('# %s::: %s\n'%(k,v))
		of1.write('# validated::: 1\n') # indicates that these results have been user-checked

		of1.write('#-\n')

		of1.write('## mz, rt, score\n')
		# get data from rows of accepted hit table
		for r in xrange(self.HT_RV_accepted_list.rowCount()):
			hitNumber = int(str(self.HT_RV_accepted_list.item(r, 0).text()))
			rt = str(self.HT_RV_accepted_list.item(r, 1).text())
			mz = str(self.HT_RV_accepted_list.item(r, 2).text())
			score = str(self.HT_RV_accepted_list.item(r, 3).text())
			of1.write('%s, %s, %s\n' %(mz, rt, score))
		of1.close()
		return

	def getXRange(self, hit, window, rangeType = 'pc'):
		'''
		Get the xrange values at +/- Xpc of a target value
		'''
		hit = float(hit)
		window = float(window)

		if rangeType == 'pc':
			xmin = hit - hit * window/100
			xmax = hit + hit * window/100
			return (xmin, xmax)

		elif rangeType == 'fixed':
			center = hit + self.mzDelta/2
			xmin = center - window
			xmax = center + window
			return (xmin, xmax)

	def getYRange(self, xData, yData, xRange, action):

		xmin, xmax = xRange

		# find indices of xData entries within xrange
		ind = np.where((xData > xmin) & (xData < xmax))

		# find y values corresponding to these points
		yDataInXR = yData[ind[0]]

		return float(np.max(yDataInXR))

	def reset(self):
		# clear contents of tablewidgets
		# clear hit list
		for i in reversed(range(self.HT_RV_hitlist.rowCount())):
			self.HT_RV_hitlist.removeRow(i)

		# clear accepted list
		for i in reversed(range(self.HT_RV_accepted_list.rowCount())):
			self.HT_RV_accepted_list.removeRow(i)

		# clear plots)
		self.HT_RV_HM_widget.clear()
		self.HT_RV_MS.clear()
		self.HT_RV_EIC.clear()
		self.makeRVPlots()
                return

        def makeRVPlots(self):

            # HM
            self.HT_RV_HM = pg.ScatterPlotItem()
            # nb labels/gridlines are applied to the widget canvas
            # rather than the scatterplotitem
            self.HT_RV_HM_widget.showGrid(x = True, y = True)
            self.HT_RV_HM_widget.showLabel('bottom', show = True)
            self.HT_RV_HM_widget.showLabel('top', show = True)
            self.HT_RV_HM_widget.showLabel('left', show = True)
            self.HT_RV_HM_widget.showLabel('right', show = True)
            self.HT_RV_HM_widget.setLabel(axis = 'bottom', text = 'm/z')
            self.HT_RV_HM_widget.setLabel(axis = 'left', text = 'Retention Time (Min)')
            self.HT_RV_HM_widget.setLabel(axis = 'top', text = '')
            self.HT_RV_HM_widget.setLabel(axis = 'right', text = '')
            self.HT_RV_HM_widget.addItem(self.HT_RV_HM)

            # EIC
            self.HT_RV_EIC.showGrid(x = True, y = True)
            self.HT_RV_EIC.showLabel('bottom', show = True)
            self.HT_RV_EIC.showLabel('left', show = True)
            self.HT_RV_EIC.showLabel('right', show = True)
            self.HT_RV_EIC.setLabel(axis = 'bottom', text = 'Retention Time (Min)')
            self.HT_RV_EIC.setLabel(axis = 'left', text = 'Intensity')
            self.HT_RV_EIC.setLabel(axis = 'top', text = '')
            self.HT_RV_EIC.setLabel(axis = 'right', text = '')

            # MS
            self.HT_RV_MS.showGrid(x = True, y = True)
            self.HT_RV_MS.showLabel('bottom', show = True)
            self.HT_RV_MS.showLabel('left', show = True)
            self.HT_RV_MS.showLabel('right', show = True)
            self.HT_RV_MS.setLabel(axis = 'bottom', text = 'm/z')
            self.HT_RV_MS.setLabel(axis = 'left', text = 'Intensity')
            self.HT_RV_MS.setLabel(axis = 'top', text = '')
            self.HT_RV_MS.setLabel(axis = 'right', text = '')

	    self.HT_RV_HM.sigClicked.connect(self.pointClicked)
            return

