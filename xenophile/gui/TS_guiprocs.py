import sys, os

import T_search
from xenophile.common import *
from xenophile.libs.targeted_PTM import results_file_parser as RFP
from xenophile.libs.targeted_PTM import targeted_PTM as TPTM
import time
import Queue
import numpy as np
import pyqtgraph as pg
from PyQt4 import QtCore, QtGui
from pyqtgraph.Point import Point

'''
Targeted search
'''
class T_search(QtGui.QDialog, T_search.Ui_Dialog):
	'''
	Targeted search methods and results viewer
	-------------------------------------------------------

	Notes:
		1) for RV tab: indices of hit objects are 1 INDEXED while rows in TableWidgets are 0 INDEXED - be careful here
	'''

	def __init__(self, darkMode = True, parent = None, runner = None):

		# set PG config options
		if not darkMode:
			pg.setConfigOption('background','w')

		super(T_search, self).__init__(parent)
		self.setupUi(self)
		self.setWindowTitle('Targeted Protein-CRM Adduct Search')

		''' create connections '''
		self.TS_gui_connections()
		self.RV_gui_connections()

		if runner:
			self.runner, self.runner_thread, self.q = runner

		''' class level lists '''
		self.mascotInputFile = None
		self.mzmlInputFile = None
		self.mzml2InputFIle = None
		self.hitimeInputFile = []
		self.GUI_input_files = [] # stores file I/O data for TS tab
		self.all_mods = [] # all variable mods detected in mascot file
		self.CRM_mods = [] # mods specified by the user as arising from CRMs
		self.HT_data = [] # ht hit objects returned from RFP - used in HM plotting

		''' set default file broser path and timer '''
		self.fileOpenDialogPath = os.path.expanduser('~')

		self.timer = QtCore.QTimer()
		self.timer.timeout.connect(self.check_response)

		''' disable mzML input until plot EIC checked '''
		self.toggle_mzml()

		''' set plot config options '''
		# HM
		self.RV_HM.showGrid(x = True, y = True)
		self.RV_HM.showLabel('bottom', show = True)
		self.RV_HM.showLabel('top', show = True)
		self.RV_HM.showLabel('left', show = True)
		self.RV_HM.showLabel('right', show = True)
		self.RV_HM.setLabel(axis = 'bottom', text = 'm/z')
		self.RV_HM.setLabel(axis = 'left', text = 'Retention Time (s)')
		self.RV_HM.setLabel(axis = 'top', text = '')
		self.RV_HM.setLabel(axis = 'right', text = '')

		# EIC
		self.RV_EIC.showGrid(x = True, y = True)
		self.RV_EIC.showLabel('bottom', show = True)
		self.RV_EIC.showLabel('left', show = True)
		self.RV_EIC.showLabel('right', show = True)
		self.RV_EIC.setLabel(axis = 'bottom', text = 'Retention Time (s)')
		self.RV_EIC.setLabel(axis = 'left', text = 'Intensity')
		self.RV_EIC.setLabel(axis = 'top', text = '')
		self.RV_EIC.setLabel(axis = 'right', text = '')

		# MS
		self.RV_spectrum.showGrid(x = True, y = True)
		self.RV_spectrum.showLabel('bottom', show = True)
		self.RV_spectrum.showLabel('left', show = True)
		self.RV_spectrum.showLabel('right', show = True)
		self.RV_spectrum.setLabel(axis = 'bottom', text = 'm/z')
		self.RV_spectrum.setLabel(axis = 'left', text = 'Intensity')
		self.RV_spectrum.setLabel(axis = 'top', text = '')
		self.RV_spectrum.setLabel(axis = 'right', text = '')

		''' set tablewidget view - remove row counter column '''
		self.RV_peptides.verticalHeader().setVisible(False)
		self.RV_accepted_peptides.verticalHeader().setVisible(False)

		''' set tablewidget selection behaviour '''
		self.T_all_mods.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
		self.T_CRM_mods.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
		self.RV_peptides.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
		self.RV_accepted_peptides.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)

		''' resize column/row contents for peptide details tableWidget '''
		self.RV_stats.resizeRowsToContents()
		self.RV_stats.resizeColumnsToContents()
		self.RV_peptides.resizeRowsToContents()
		self.RV_peptides.resizeColumnsToContents()

		# stretch last column to fill widget space - need a better way of doing shis :(
		self.RV_stats.horizontalHeader().setStretchLastSection(True)
		self.RV_peptides.horizontalHeader().setStretchLastSection(True)
		self.T_all_mods.horizontalHeader().setStretchLastSection(True)
		self.T_CRM_mods.horizontalHeader().setStretchLastSection(True)
		self.T_HT_files.horizontalHeader().setStretchLastSection(True)
		self.RV_accepted_peptides.horizontalHeader().setStretchLastSection(True)

		# testing
		self.HT_score_cutoff.setText('0')
		self.T_HT_rtGap.setText('30')
		self.T_HT_mzGap.setText('6')
		self.T_NMD.setText('6.0201')
		self.T_HT_peptide_correlation_RT.setText('2')
		self.T_HT_peptide_correlation_mz.setText('0.5')
		self.T_eicWidth.setText('0.03')

	def TS_gui_connections(self):
		''' create connections for TPS search tab '''
		self.T_HT_add.clicked.connect(lambda: self.TS_select_file('TS_HT_entry_button'))
		self.T_browse_mascot.clicked.connect(lambda: self.TS_select_file('TS_mascot_entry_button'))
		self.T_browse_mzML.clicked.connect(lambda: self.TS_select_file('TS_mzml_entry_button'))
		self.T_browse_mzML2.clicked.connect(lambda: self.TS_select_file('TS_mzml_entry_button2'))
		self.T_HT_remove.clicked.connect(self.remove_HT_file_from_table)
		self.T_mods_add.clicked.connect(self.add_highlighted_mod_to_CRM_list)
		self.T_mods_remove.clicked.connect(self.remove_highlighted_mod_from_CRM_list)
		self.T_run.clicked.connect(self.run_targeted_search)
		self.T_outputFile_button.clicked.connect(self.set_output_file_dir)
		self.T_plot_EICs.toggled.connect(self.toggle_mzml)
		return

	def check_response(self):
		'''
		On timer tick - check response from NTPS
			- if process returns done, load results into RV GUI tab
		'''
		while not self.q.empty():
			update = self.q.get()
			if update == 'done':
				time.sleep(0.1)
				self.runner.p.terminate()
				self.timer.stop()
				self.runner_thread.exit()
				print('\nResponse returned to guiprocs')
			else:
				print update
		return

	def set_output_file_dir(self):
		outputFileName = str(QtGui.QFileDialog.getSaveFileName(self, 'Set Output File', self.fileOpenDialogPath))
		self.fileOpenDialogPath = os.path.dirname(outputFileName)
		self.T_outputFile_text.setText(outputFileName)
		return

	def toggle_mzml(self):
		state = self.T_plot_EICs.isChecked()
		self.T_browse_mzML.setEnabled(state)
		self.T_mzM_file.setEnabled(state)
		self.T_mzML_file2.setEnabled(state)
		self.T_browse_mzML2.setEnabled(state)
		return

	def TS_select_file(self, caller):
		''' open file browser and init file handling actions '''
		if caller == 'TS_HT_entry_button':
			TS_HT_file = str(QtGui.QFileDialog.getOpenFileName(self, 'Select HiTIME File', self.fileOpenDialogPath))
			self.fileOpenDialogPath = os.path.dirname(str(TS_HT_file))
			self.hitimeInputFile.append(TS_HT_file)

			# add HT file to TableWidget
			# find number of rows
			rowCount = self.T_HT_files.rowCount()

			# add a new row to the end of the table
			self.T_HT_files.insertRow(rowCount)

			# insert data
			self.T_HT_files.setItem(rowCount, 0, QtGui.QTableWidgetItem(str(TS_HT_file)))
			self.T_HT_files.setItem(rowCount, 1, QtGui.QTableWidgetItem(str('<charge>')))

			# centre text for charge - leave name as left justified
			self.T_HT_files.item(rowCount,1).setTextAlignment(QtCore.Qt.AlignCenter)

			# resize column contents
			self.T_HT_files.resizeRowsToContents()
			self.T_HT_files.resizeColumnsToContents()

			table_header = self.T_HT_files.horizontalHeader()
			table_header.setStretchLastSection(True)

		if caller == 'TS_mascot_entry_button':
			TS_mascot_file = str(QtGui.QFileDialog.getOpenFileName(self, 'Select Mascot File', self.fileOpenDialogPath))
			self.fileOpenDialogPath = os.path.dirname(str(TS_mascot_file))
			self.mascotInputFile = TS_mascot_file
			self.T_mascot_file.setText(TS_mascot_file)
			self.get_mascot_mods(TS_mascot_file)

		if caller == 'TS_mzml_entry_button':
			TS_mzml_file = str(QtGui.QFileDialog.getOpenFileName(self, 'Select MS1 mzML File', self.fileOpenDialogPath))
			self.fileOpenDialogPath = os.path.dirname(str(TS_mzml_file))
			self.T_mzM_file.setText(TS_mzml_file)
		if caller == 'TS_mzml_entry_button2':
			TS_mzml_file2 = str(QtGui.QFileDialog.getOpenFileName(self, 'Select MS2 mzML File', self.fileOpenDialogPath))
			self.fileOpenDialogPath = os.path.dirname(str(TS_mzml_file2))
			self.T_mzML_file2.setText(TS_mzml_file2)
		return

	def remove_HT_file_from_table(self):
		highlighted_row = self.T_HT_files.selectionModel().selectedRows()[0].row()
		self.T_HT_files.removeRow(highlighted_row)
		return

	def get_mascot_mods(self, ifpath):
		''' retrieve PTMs from mascot output file and add to tablewidget '''

		# Get var mods ---> returns list of varMod objects.
		# 	Attributes are: modIndex, modName, modDelta, modNeutralLoss
		var_mods = get_mods_from_mascot_file(str(ifpath))

		# add var_mods to class level list
		for i in var_mods: self.all_mods.append(i)

		for mod in var_mods:
			# add row to all mods tablewidget
			rowCount = self.T_all_mods.rowCount()
			self.T_all_mods.insertRow(rowCount)

			# insert data
			self.T_all_mods.setItem(rowCount, 0, QtGui.QTableWidgetItem(str(mod.modName)))
			self.T_all_mods.setItem(rowCount, 1, QtGui.QTableWidgetItem(str(mod.modDelta)))

			# centre text for mass - leave name as left justified
			self.T_all_mods.item(rowCount,1).setTextAlignment(QtCore.Qt.AlignCenter)

			# resize column contents
			self.T_all_mods.resizeRowsToContents()
			self.T_all_mods.resizeColumnsToContents()
		return

	def add_highlighted_mod_to_CRM_list(self):
		''' add modification highlighted in all mods table to CRM mods table '''
		# get highlighted row
		highlighted_row = self.T_all_mods.selectionModel().selectedRows()[0].row()

		# get mod from all mods class level list
		mod = self.all_mods[highlighted_row]

		# determine number of rows CRM table
		rowCount = self.T_CRM_mods.rowCount()

		# add a new row
		self.T_CRM_mods.insertRow(rowCount)

		# insert data
		self.T_CRM_mods.setItem(rowCount, 0, QtGui.QTableWidgetItem(str(mod.modName)))
		self.T_CRM_mods.setItem(rowCount, 1, QtGui.QTableWidgetItem(str(mod.modDelta)))

		# centre text for mass - leave name as left justified
		self.T_CRM_mods.item(rowCount,1).setTextAlignment(QtCore.Qt.AlignCenter)

		# resize column contents
		self.T_CRM_mods.resizeRowsToContents()
		self.T_CRM_mods.resizeColumnsToContents()
		self.T_CRM_mods.horizontalHeader().setStretchLastSection(True)
		return

	def remove_highlighted_mod_from_CRM_list(self):
		''' remove highlighted row from CRM mods table '''
		# get highlighted row
		highlighted_row = self.T_CRM_mods.selectionModel().selectedRows()[0].row()
		# remove row
		self.T_CRM_mods.removeRow(highlighted_row)
		return

	def retrieve_mod_index(self):
		'''
		Get mascot index associated with the entries in the
		mascot CRM target list
		'''
		rowCount = self.T_CRM_mods.rowCount()

		for r in xrange(rowCount):
			modName = str(self.T_CRM_mods.item(r, 0).text())

			# get corresponding entry from varMods list
			for vm in self.all_mods:
				if modName == vm.modName:
					self.CRM_mods.append(vm)
		return

	def run_targeted_search(self):
		print 'ts guiprocs run targeted: %s' % QtCore.QThread.currentThreadId()
		# pull options from GUI
		mascotInput = str(self.T_mascot_file.text())

		htInput = []
		for row in range(int(self.T_HT_files.rowCount())):
			inFile = str(self.T_HT_files.item(row, 0).text())
			charge = str(self.T_HT_files.item(row, 1).text())
			htInput.append([inFile, charge])

		scoreCutoff = float(self.HT_score_cutoff.text())
		rtGap = float(self.T_HT_rtGap.text())
		mzGap = float(self.T_HT_mzGap.text())
		getIndirect = self.T_get_indirect.isChecked()
		outputfile = str(self.T_outputFile_text.text())
		targetMods = self.retrieve_mod_index()
		modIndices = [i.modIndex for i in self.CRM_mods]
		plotData = self.T_plot_EICs.isChecked()
		neutralMD = float(str(self.T_NMD.text()))
		correlationRT = float(self.T_HT_peptide_correlation_RT.text())
		correlationMZ = float(self.T_HT_peptide_correlation_mz.text())

		if plotData:
			mzmlInput = str(self.T_mzM_file.text())
			mzmlMS2File = str(self.T_mzML_file2.text())
			eicWidth = float(self.T_eicWidth.text())
		else:
			mzmlInput = None
			mzmlMS2File = None
			eicWidth = None

		# produce argument dictionary
		args = {
			'outputfile' : outputfile,
			'htInput' : htInput,
			'mascotInput' : mascotInput,
			'getIndirect' : getIndirect,
			'scoreCutoff' : scoreCutoff,
			'rtGap' : rtGap,
			'mzGap' : mzGap,
			'modIndices' : modIndices,
			'mzmlFile' : mzmlInput,
			'mzmlMS2File' : mzmlMS2File,
			'plotResults' : plotData,
			'mzDelta' : neutralMD,
			'eicWidth' : eicWidth,
			'correlationRT': correlationRT,
			'correlationMZ' : correlationMZ
			#'htCharge' : htCharge,
		}

		run_job(self.runner, self.runner_thread, self.q, TPTM.params, args)
		self.timer.start(500)

	'''
	RV tab functions
	'''
	def RV_init_routines(self):
		# add empty rows to RV_stats tableview widget
		rowPosition = self.RV_stats.rowCount()
		for i in range(10):
			self.RV_stats.insertRow(rowPosition)

		# add data labels to row 1
		self.RV_stats.setItem(0, 0, QtGui.QTableWidgetItem("text1"))

		return

	def RV_gui_connections(self):
		self.RV_load_results.clicked.connect(self.parse_results_file)
		self.RV_peptides.itemSelectionChanged.connect(lambda: self.update_RV_plots(None))
		self.RV_up.clicked.connect(lambda: self.update_RV_plots('up'))
		self.RV_down.clicked.connect(lambda: self.update_RV_plots('down'))
		self.RV_MS_RB.toggled.connect(lambda: self.update_RV_plots(None))
		self.RV_MSMS_RB.toggled.connect(lambda: self.update_RV_plots(None))
		self.RV_accept.clicked.connect(self.add_peptide_to_accepted_list)
		self.RV_reject.clicked.connect(lambda: self.update_RV_plots('down'))
		self.RV_accept.clicked.connect(lambda: self.update_RV_plots('down'))
		self.RV_remove.clicked.connect(self.remove_peptide_from_accepted_list)
		self.RV_save.clicked.connect(self.write_output_file)
		return

	def parse_results_file(self):
		''' called when results file is loaded '''

		print('parsing results file')
		# need to spawn file selection panel to retrieve results files
		# parse results file
		inFile = QtGui.QFileDialog.getOpenFileName(self, 'Select Targeted Search Output File', self.fileOpenDialogPath)
		self.fileOpenDialogPath = os.path.dirname(str(inFile))
		self.mascot_hits, self.HT_data = RFP.run(ifile = inFile) # use testing file for the moment

		print('finished parsing results file')
		# set num rows in tableview widget
		#self.RV_peptides.setRowCount(len(self.mascot_hits))

		# populate table
		for i, hit in enumerate(self.mascot_hits):
			self.RV_peptides.insertRow(i)
			self.RV_peptides.setItem(i, 0, QtGui.QTableWidgetItem(str(hit.index)))
			self.RV_peptides.setItem(i, 1, QtGui.QTableWidgetItem(hit.Type))
			self.RV_peptides.setItem(i, 2, QtGui.QTableWidgetItem('%.2f' %float(hit.HT_rt)))
			self.RV_peptides.setItem(i, 3, QtGui.QTableWidgetItem('%.2f' %float(hit.mz)))

			# center text
			self.RV_peptides.item(i,0).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_peptides.item(i,1).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_peptides.item(i,2).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_peptides.item(i,3).setTextAlignment(QtCore.Qt.AlignCenter)

		# resize column contents
		self.RV_peptides.resizeRowsToContents()
		self.RV_peptides.resizeColumnsToContents()

		# stretch last column to fill widget space
		self.RV_peptides.horizontalHeader().setStretchLastSection(True)
		self.RV_stats.horizontalHeader().setStretchLastSection(True)

		# set first row as highlighted
		self.RV_peptides.selectRow(0)

		# plot heatmap data - NB: do this once then adjust view as needed
		self.plot_TS_RV_HM()

		# add data plots for first row
		self.update_RV_plots(None)

		return

	def plot_TS_RV_HM(self):
		'''
		Plot heatmap data derived from results file.
		NB - heatmap data is not explicitely written to the output
			--> points associated with peptide hits can be derived from the hit data but the full heatmap data is not necessarily present
			--> could either append HM data to the results file or ask for HM file EIC_specification
				---> list of HM peak maxima is probably small so append to the results file
		'''

		self.RV_HM.plot(x = np.asarray([x.mz for x in self.HT_data]),
						y = np.asarray([x.rt for x in self.HT_data]),
						symbol = 'o',
						pen = None)
		return

	def update_RV_plots(self, direction):
		'''
		Update data analysis plots and stats when
			1) up or down button is clicked, or
			2) a new peptide entry is clicked in the results list
			3) MS spectrum type is changed using radio buttons

		Parameters:
			direction:
				'up' if up button clicked
				'down' if down dubbon clicked
				'None' if new entry clicked in results list

		'''

		#print('\n\n')
		# Set default radio button status if None
		if self.RV_MS_RB.isChecked() == False and self.RV_MSMS_RB.isChecked() == False:
			self.RV_MS_RB.setChecked(True)

		# get selected row from RV_peptides table
		if direction is not None:
			# direction != None when up or down buttons are clicked
			index = self.RV_peptides.selectionModel().selectedRows()[0].row()

			if direction == 'up':
				self.RV_peptides.selectRow(index-1)
			elif direction == 'down':
				self.RV_peptides.selectRow(index+1)

		# ensure that a row is highlighted
		iaARowHighlighted = self.RV_peptides.selectionModel().selectedRows()
		if len(iaARowHighlighted) == 0: return

		# fomd index of highlighted row
		highlighted_row = self.RV_peptides.selectionModel().selectedRows()[0].row()

		# get value in first column corresponding to the hit object attribute index
		item = self.RV_peptides.item(highlighted_row, 0).text()
		#print('Item at row %s, column %s is: %s' %(highlighted_row, 0, item))
		#print('selected indices are: %s' %([x.row() for x in self.RV_peptides.selectionModel().selectedRows()]))
		#aprint('highlighted row is: %s' %highlighted_row)

		# get matching hit
		hit = None
		for x in self.mascot_hits:
			if str(x.index) == str(item):
				#print(str(highlighted_row), str(x.index))
				hit = x
				break

		if hit is not None:
			# Update peptide details table
			self.RV_stats.setItem(0, 0, QtGui.QTableWidgetItem(hit.sequence))
			self.RV_stats.setItem(0, 1, QtGui.QTableWidgetItem(hit.charge))
			self.RV_stats.setItem(0, 2, QtGui.QTableWidgetItem(hit.MC))
			self.RV_stats.setItem(0, 3, QtGui.QTableWidgetItem(hit.mod_string))
			self.RV_stats.setItem(0, 4, QtGui.QTableWidgetItem(hit.pep_score)) # score
			self.RV_stats.setItem(0, 5, QtGui.QTableWidgetItem(hit.pep_identity_score)) # identity
			self.RV_stats.setItem(0, 6, QtGui.QTableWidgetItem(hit.pep_homology_score)) # homology
			self.RV_stats.setItem(0, 7, QtGui.QTableWidgetItem(hit.prot_matches))
			self.RV_stats.setItem(0, 8, QtGui.QTableWidgetItem(hit.prot_acc))
			self.RV_stats.setItem(0, 9, QtGui.QTableWidgetItem(hit.mz))
			self.RV_stats.setItem(0, 10, QtGui.QTableWidgetItem(hit.rt))
			self.RV_stats.setItem(0, 11, QtGui.QTableWidgetItem(hit.HT_mz))
			self.RV_stats.setItem(0, 12, QtGui.QTableWidgetItem(hit.HT_rt))

			# center text - this could be neater...
			self.RV_stats.item(0,0).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,1).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,2).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,3).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,4).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,5).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,6).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,7).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,8).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,9).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,10).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,11).setTextAlignment(QtCore.Qt.AlignCenter)
			self.RV_stats.item(0,12).setTextAlignment(QtCore.Qt.AlignCenter)

			# resize columns
			self.RV_stats.resizeRowsToContents()
			self.RV_stats.resizeColumnsToContents()

			# stretch last column to fill widget space
			self.RV_stats.horizontalHeader().setStretchLastSection(True)

			# Prepare data for plotting
			# EIC data
			RT_data = np.asarray(self.clean_data(hit.EIC_RT))
			Ilight_data = np.asarray(self.clean_data(hit.EIC_light))
			Iheavy_data = np.asarray(self.clean_data(hit.EIC_heavy))

			# MS data
			# peptide
			pep_ms1_mz = np.asarray(self.clean_data(hit.pep_MS_mz), dtype = 'float32') + 0.05
			pep_ms1_int = np.asarray(self.clean_data(hit.pep_MS_int))

			# print pep_ms1_mz.shape
			# print pep_ms1_int.shape

			# HT hit
			ht_ms1_mz = np.asarray(self.clean_data(hit.HT_MS_mz))
			ht_ms1_int = np.asarray(self.clean_data(hit.HT_MS_int))

			# print ht_ms1_mz.shape
			# print ht_ms1_int.shape

			# MS/MS data
			ms2_mz = np.asarray(self.clean_data(hit.pep_MS2_mz))
			ms2_int = np.asarray(self.clean_data(hit.pep_MS2_int))

			# clear current plots
			self.RV_EIC.clear()
			self.RV_spectrum.clear()

			# remove old plot legends and add new ones
			try:
				self.eicLGD.scene().removeItem(self.eicLGD)
				self.msLGD.scene().removeItem(self.msLGD)
			except:
				pass

			self.eicLGD = self.RV_EIC.addLegend()
			self.msLGD = self.RV_spectrum.addLegend()

			# self.RV_HM.clear()

			# EIC
			self.RV_EIC.plot(x = RT_data, y = Ilight_data, pen = 'b', name = 'Light')
			self.RV_EIC.plot(x = RT_data, y = Iheavy_data, pen = 'r', name = 'Heavy')

			# MS spectrum
			if self.RV_MS_RB.isChecked():

				pep_ms1_mz, pep_ms1_int = zero_fill(pep_ms1_mz, pep_ms1_int)
				self.RV_spectrum.plot(x = pep_ms1_mz, y = pep_ms1_int, pen = 'b', name = 'Pep, RT = %s, m/z = %s' %(hit.rt, hit.mz))

				ht_ms1_mz, ht_ms1_int = zero_fill(ht_ms1_mz, ht_ms1_int)
				self.RV_spectrum.plot(x = ht_ms1_mz, y = ht_ms1_int, pen = 'r', name = 'HT, RT = %s, m/z = %s' %(hit.HT_rt, hit.HT_mz))
				xRangeMS = self.getXRange(hit.HT_mz, 10, rangeType = 'fixed')
				yRangeMS = self.getYRange(ht_ms1_mz, ht_ms1_int, xRangeMS)

			elif self.RV_MSMS_RB.isChecked():
				xRangeMS = (float(min(ms2_mz)), float(max(ms2_mz)))
				yRangeMS = (float(min(ms2_int)), float(max(ms2_int)))
				self.RV_spectrum.plot(x=np.repeat(ms2_mz, 3), y=np.dstack((np.zeros(ms2_int.shape[0]), ms2_int, np.zeros(ms2_int.shape[0]))).flatten(), connect='all', pen = 'b')
				self.RV_spectrum.setLimits(xMin = xRangeMS[0], xMax = xRangeMS[1], yMin = yRangeMS[0], yMax = yRangeMS[1])
				self.RV_spectrum.setRange(xRange = xRangeMS, yRange = yRangeMS)


			# set ranges
			xRangeRT = self.getXRange(hit.HT_rt, 10, rangeType = 'pc')
			yRangeRTL = self.getYRange(RT_data, Ilight_data, xRangeRT)
			yRangeRTH = self.getYRange(RT_data, Iheavy_data, xRangeRT)
			yRangeRT = max([yRangeRTL, yRangeRTH])

			yRangeRT = ( 0 - yRangeRT * 0.05, yRangeRT)
			self.RV_EIC.setLimits(xMin = float(min(RT_data)), xMax = float(max(RT_data)), yMin = yRangeRT[0], yMax = yRangeRT[1]*1.5)
			self.RV_EIC.setRange(xRange = xRangeRT, yRange = yRangeRT)

			if not self.RV_MSMS_RB.isChecked():
				yRangeMS = ( 0 - yRangeMS * 0.05, yRangeMS)
				self.RV_spectrum.setLimits(xMin = float(min(ht_ms1_mz)), xMax = float(max(ht_ms1_mz)), yMin = yRangeRT[0], yMax = yRangeMS[1]*1.5)
				self.RV_spectrum.setRange(xRange = xRangeMS, yRange = yRangeMS)


		return

	def clean_data(self, data):
		clean_data = []
		for i in data:
			try:
				clean_data.append(float(i))
			except:
				#print(i)
				pass

		return clean_data

	def add_peptide_to_accepted_list(self):
		''' add selected peptide to accepted hits list '''
		# get index of highlighted row
		highlighted_row = self.RV_peptides.selectionModel().selectedRows()[0].row()
		item = self.RV_peptides.item(highlighted_row, 0).text()
		# get matching hit
		hit = None
		for x in self.mascot_hits:
			if str(x.index) == str(item):
				#print(str(highlighted_row), str(x.index))
				hit = x
				break

		# add peptie to accepted list
		rowCount = self.RV_accepted_peptides.rowCount()
		self.RV_accepted_peptides.insertRow(rowCount)

		# add data
		self.RV_accepted_peptides.setItem(rowCount, 0, QtGui.QTableWidgetItem(str(hit.index)))
		self.RV_accepted_peptides.setItem(rowCount, 1, QtGui.QTableWidgetItem(hit.sequence))
		self.RV_accepted_peptides.setItem(rowCount, 2, QtGui.QTableWidgetItem(hit.prot_acc))

		# center text
		self.RV_accepted_peptides.item(rowCount,0).setTextAlignment(QtCore.Qt.AlignCenter)
		self.RV_accepted_peptides.item(rowCount,1).setTextAlignment(QtCore.Qt.AlignCenter)
		self.RV_accepted_peptides.item(rowCount,2).setTextAlignment(QtCore.Qt.AlignCenter)

		# resize columns
		self.RV_accepted_peptides.resizeRowsToContents()
		self.RV_accepted_peptides.resizeColumnsToContents()

		# stretch last column to fill widget space
		table_header = self.RV_accepted_peptides.horizontalHeader()
		table_header.setStretchLastSection(True)
		return

	def remove_peptide_from_accepted_list(self):
		''' remove selected peptide to accepted hits list '''
		# get highlighted row
		highlighted_row = self.RV_accepted_peptides.selectionModel().selectedRows()[0].row()
		# remove row
		self.RV_accepted_peptides.removeRow(highlighted_row)
		return

	def write_output_file(self):
		outfile = str(QtGui.QFileDialog.getSaveFileName(self, 'Select Output File', self.fileOpenDialogPath))
		of1 = open(outfile, 'wt')

		rc = self.RV_accepted_peptides.rowCount()

		for i in range(rc):
			result = None
			rowIndex = str(self.RV_accepted_peptides.item(i, 0).text())
			for j in self.mascot_hits:
				if str(rowIndex) == str(j.index):
					result = j
					break

			if not result: continue

			print >> of1, 'Results for Hit {}'.format(result.index)
			print >> of1, '..................'.format()
			print >> of1, 'PROTEIN_ACC = {}'.format(result.prot_acc)
			print >> of1, 'PROTEIN_INDEX = {}'.format(result.prot_index)
			print >> of1, 'PROTEIN_MATCHES = {}'.format(result.prot_matches)
			print >> of1, 'PEPTIDE_MZ = {}'.format(result.mz)
			print >> of1, 'PEPTIDE_RT = {}'.format(result.rt)
			print >> of1, 'HITIME_MZ = {}'.format(result.HT_mz)
			print >> of1, 'HITIME_RT = {}'.format(result.HT_rt)
			print >> of1, 'PEPTIDE_CHARGE = {}'.format(result.charge)
			print >> of1, 'PEPTIDE_MISSED_CLEAVAGE = {}'.format(result.MC)
			print >> of1, 'PEPTIDE_SEQUENCE = {}'.format(result.sequence)
			print >> of1, 'PEPTIDE_VARMODS = {}'.format(result.mod_string)
			print >> of1, 'PEPTIDE_SCORE = {}'.format(result.pep_score)
			print >> of1, 'PEPTIDE_IDENTITY_SCORE = {}'.format(result.pep_identity_score)
			print >> of1, 'PEPTIDE_HOMOLOGY_SCORE = {}'.format(result.pep_homology_score)
			print >> of1, ''.format()
			print >> of1, ''.format()

		of1.close()
		return

	def getXRange(self, hit, window, rangeType = 'pc'):
		'''
		Get the xrange values at +/- Xpc of a target value
		'''
		hit = float(hit)
		window = float(window)
		self.mzDelta = 6.0201 # TODO, this can't be hard-coded
		if rangeType == 'pc':
			xmin = hit - hit * window/100
			xmax = hit + hit * window/100
			return (xmin, xmax)

		elif rangeType == 'fixed':
			center = hit + self.mzDelta/2
			xmin = center - window
			xmax = center + window
			return (xmin, xmax)

	def getYRange(self, xData, yData, xRange):

		xmin, xmax = xRange

		# find indices of xData entries within xrange
		ind = np.where((xData > xmin) & (xData < xmax))

		# find y values corresponding to these points
		yDataInXR = yData[ind]

		# return max Y value
		return float(np.max(yDataInXR))
