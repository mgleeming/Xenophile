
import sys, os

# Add JSME directory to PATH
sys.path.append(
	os.path.normpath(
		os.path.join(
			os.path.dirname(__file__), os.pardir,'libs','JSME'
			)
		)
	)

#print os.listdir(os.path.join(os.path.dirname(__file__), os.pardir,'libs','JSME'))

from xenophile.common import *
import NT_search
import xenophile.libs.non_targeted_PTM.non_targeted_PTM_identifier as NTPI

from PyQt4 import QtCore, QtGui
from PyQt4.QtWebKit import QWebView # for JSME
import pyqtgraph as pg
from pyqtgraph.Point import Point

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from pyteomics import mgf, pepxml, mass
from pyteomics import mgf, pepxml, mass
from multiprocessing import Process, Queue

'''
NTPS search tab
'''
class NT_search(QtGui.QDialog, NT_search.Ui_Dialog):
	'''
	Non-targeted search methods and results viewer
	'''

	'''
	NOTES:

	It seems that the m/z value determined by selecting local maxima of HT peaks is slightly inaccurate
		- this introduces errors into the calculation of deltaMZ which is taken as the CRM mass
		- As a workaround - find the MGF entry associated with the HT hit and use this m/z value which should be more accurate

	For NTPS to work - HT hits must have an associated MS2 spectrum

	'''

	def __init__(self, darkMode = True, parent = None, runner = None):

		# set PG config options
		if not darkMode:
			pg.setConfigOption('background','w')

		self.darkMode = darkMode

		super(NT_search, self).__init__(parent)
		self.setupUi(self)

		self.setWindowTitle('Non-targeted Reactive Metabolite Search')

		if runner:
			self.runner, self.runner_thread, self.q = runner

		''' create connections '''
		self.NTPS_gui_connections()
		self.NTPS_RV_gui_connections()

		''' set temp values for debugging '''
		self.NTPS_smiles_entry.setText('CC(=O)Nc1ccc(O)cc1')
		self.NTPS_min_HT_score.setText('15')

		self.NTPS_HT_MS2_mz_offset.setText('0.5')
		self.NTPS_HT_MS2_rt_offset.setText('2')


		''' set timer '''
		self.timer = QtCore.QTimer()
		self.timer.timeout.connect(self.check_response)

		''' init multiprocessing queue'''
		self.q = Queue(maxsize = 0)

		''' multiprocess job '''
		self.job = None

		''' class level data '''
		self.molecule = [] # store RDkit descriptions of mol and fragments
		self.structures = [] # filenames of fragment structures
		self.GUI_input_files = [] # stores file I/O data for HTS tab
		self.mascot_hits = [] # mascot hits imported into RV tab
		self.hits = [] # NTPS search results returned and used in RV tab
		self.heatmap = [] # HT local maxima returned by NTPS - list of HT objects

		self.fileOpenDialogPath = os.path.expanduser('~')

		''' set tablewidget view - remove row counter column '''
		self.NTPS_RV_hits.verticalHeader().setVisible(False)
		self.NTPS_RV_formulae.verticalHeader().setVisible(False)
		self.NTPS_RV_accepted_hits.verticalHeader().setVisible(False)

		''' set plot config options '''
		# HM
		self.NTPS_RV_HM.showGrid(x = True, y = True)
		self.NTPS_RV_HM.showLabel('bottom', show = True)
		#self.NTPS_RV_HM.showLabel('top', show = True)
		self.NTPS_RV_HM.showLabel('left', show = True)
		self.NTPS_RV_HM.showLabel('right', show = True)
		self.NTPS_RV_HM.setLabel(axis = 'bottom', text = 'm/z')
		self.NTPS_RV_HM.setLabel(axis = 'left', text = 'Retention Time (s)')
		self.NTPS_RV_HM.setLabel(axis = 'top', text = '')
		self.NTPS_RV_HM.setLabel(axis = 'right', text = '')

		# MS
		self.NTPS_RV_MS.showGrid(x = True, y = True)
		self.NTPS_RV_MS.showLabel('bottom', show = True)
		self.NTPS_RV_MS.showLabel('left', show = True)
		self.NTPS_RV_MS.showLabel('right', show = True)
		self.NTPS_RV_MS.setLabel(axis = 'bottom', text = 'm/z')
		self.NTPS_RV_MS.setLabel(axis = 'left', text = 'Intensity')
		self.NTPS_RV_MS.setLabel(axis = 'top', text = '')
		self.NTPS_RV_MS.setLabel(axis = 'right', text = '')

		# QPen for drawing theoretical MS/MS peaks
		self.qp = QtGui.QPen()
		self.qp.setColor(QtGui.QColor(255,0,0,120)) # values are R, G, B, A - A = transparency > ranges from 0 (fully tranparent) to 255 (fully opaque)
		self.qp.setStyle(QtCore.Qt.DashLine)

		''' set TableWidget selection models '''	## NEED TO LIMIT SELECTION TO ONE ROW AT A TIME
		self.NTPS_RV_hits.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
		self.NTPS_RV_accepted_hits.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
		self.NTPS_RV_formulae.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)

		# resize column contents
		self.NTPS_RV_hits.resizeRowsToContents()
		self.NTPS_RV_hits.resizeRowsToContents()
		self.NTPS_RV_accepted_hits.resizeRowsToContents()
		self.NTPS_RV_accepted_hits.resizeRowsToContents()
		self.NTPS_RV_formulae.resizeRowsToContents()
		self.NTPS_RV_formulae.resizeRowsToContents()
		self.NTPS_mod_comp_table.resizeRowsToContents()
		self.NTPS_mod_comp_table.resizeRowsToContents()

		# stretch last column to fill widget space
		self.NTPS_RV_hits.horizontalHeader().setStretchLastSection(True)
		self.NTPS_RV_accepted_hits.horizontalHeader().setStretchLastSection(True)
		self.NTPS_RV_formulae.horizontalHeader().setStretchLastSection(True)
		self.NTPS_mod_comp_table.horizontalHeader().setStretchLastSection(True)

		# cleanup old structure images
		self.cleanup()

	def set_tablewidget_column_geometries(self, widgets):
		'''
		set tablewidget row/column size properties
			- sets resize row and column
			- sets setStretchLastSection
		'''
		for item in widgets:
			item.resizeRowsToContents()
			item.resizeRowsToContents()
			table_header = item.horizontalHeader().setStretchLastSection(True)
		return

	def NTPS_gui_connections(self):
		self.NTPS_HT_entry_button.clicked.connect(lambda: self.NTPS_select_file('NTPS_HT_entry_button'))
		self.NTPS_mascot_entry_button.clicked.connect(lambda: self.NTPS_select_file('NTPS_mascot_entry_button'))
		self.NTPS_draw_smiles.clicked.connect(self.JSME_draw)
		self.NTPS_gen_frags.clicked.connect(self.generate_fragemnts)
		self.NTPS_run.clicked.connect(self.run_non_targeted_search)
		self.NTPS_frag_list.itemSelectionChanged.connect(lambda: self.update_structure_viewer(None))
		self.NTPS_frag_list_up_button.clicked.connect(lambda: self.update_structure_viewer('up'))
		self.NTPS_frag_list_down_button.clicked.connect(lambda: self.update_structure_viewer('down'))
		self.NTPS_add_element_row.clicked.connect(self.ntps_AddElementRow)
		self.NTPS_mod_comp_table.itemChanged.connect(self.ntps_updateMzBand)
		self.NTPS_output_browse.clicked.connect(self.ntps_select_output)
		return

	def launch_rangefinder(self):
		dialog = RangeFinder(darkMode = self.darkMode)
		dialog.findRanges(self)
		return

	def ntps_AddElementRow(self):
		rc = self.NTPS_mod_comp_table.rowCount()
		self.NTPS_mod_comp_table.insertRow(rc)
		return

	def NTPS_select_file(self, caller):
		''' TODO - cleanup '''
		'''
		NOTE - default file open location is the APAP lumos data directory - MUST change this before publishing
		'''
		if caller == 'NTPS_HT_entry_button':
			NTPS_HT_file = QtGui.QFileDialog.getOpenFileName(self, 'Select HiTIME File', self.fileOpenDialogPath)
			self.fileOpenDialogPath = os.path.dirname(str(NTPS_HT_file))
			#self.get_NTPS_input_file(NTPS_HT_file, 'HT')
			self.ntps_htInputFile = str(NTPS_HT_file)
			self.NTPS_HT_entry.setText(self.ntps_htInputFile)

		if caller == 'NTPS_mascot_entry_button':
			NTPS_mascot_file = QtGui.QFileDialog.getOpenFileName(self, 'Select Mascot File', self.fileOpenDialogPath)
			self.fileOpenDialogPath = os.path.dirname(str(NTPS_mascot_file))
			#self.get_NTPS_input_file(NTPS_mascot_file, 'mascot')
			self.ntps_mascotInputFile = str(NTPS_mascot_file)
			self.NTPS_mascot_entry.setText(self.ntps_mascotInputFile)
		return

	def JSME_draw(self):
		'''
		Call JSME HTML file and load into webview widget
		'''

		path_to_JSME = os.path.normpath(
							os.path.join(
								os.path.dirname(__file__), os.pardir,'libs','JSME', 'JSME_HTML_test_qt_int.html'
								)
							)

		# And a window
		self.win = QtGui.QWidget()
		self.win.setWindowTitle('Draw input molecule')

		# Create and fill a QWebView
		self.view = QWebView(self.win)
		self.view.setFixedSize(660,500)
		self.view.load(QtCore.QUrl(path_to_JSME))
		self.button = QtGui.QPushButton('Submit', parent = self.win)
		self.button.move(0,490)
		self.button.setFixedSize(660,35)
		self.button.clicked.connect(self.JSME_get_smiles)
		self.button.clicked.connect(self.win.close)
		self.win.show()
		return

	# Interact with the HTML page by calling the completeAndReturnName
	# function; print its return value to the console
	def JSME_get_smiles(self):
		'''
		Get the SMILES string of the molecule entered into JSME on button click
		'''
		self.frame = self.view.page().mainFrame()
		self.variant = self.frame.evaluateJavaScript('jsmeApplet.smiles();')
		self.a = str(self.variant.toString())
		self.NTPS_smiles_entry.setText(self.a)
		return

	def generate_fragemnts(self):
		'''
		Disconnect rotatable bond of SMILES specified molecule
		Save fragments calculate formulae
		Load fragments into structure viewer widget
		'''

		# Clear class level list
		del self.molecule[:]

		# Get text input from LineEdit widget
		self.smile = str(self.NTPS_smiles_entry.text())

		# Append cStructure mol to class level list
		self.molecule.append(cStructure('precursor', self.smile, 'precursor.png'))

		# SMARTS string for bond disconnection
		self.patt = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]')

		# Init parent mol object
		self.mol = Chem.MolFromSmiles(self.smile)

		# find the rotatable bonds
		self.bonds = self.mol.GetSubstructMatches(self.patt)

		# create an editable molecule, break the bonds, and add dummies:
		self.all_smis = [self.smile]

		# disconnect rotatable bonds
		for a,b in self.bonds:
			self.em = Chem.EditableMol(self.mol)
			self.nAts = self.mol.GetNumAtoms()
			self.em.RemoveBond(a,b)
			self.em.AddAtom(Chem.Atom(0))
			self.em.AddBond(a,self.nAts,Chem.BondType.SINGLE)
			self.em.AddAtom(Chem.Atom(0))
			self.em.AddBond(b,self.nAts+1,Chem.BondType.SINGLE)
			self.nAts+=2
			self.p = self.em.GetMol()
			Chem.SanitizeMol(self.p)
			self.smis = [Chem.MolToSmiles(x,True) for x in Chem.GetMolFrags(self.p,asMols=True)]
			for self.smi in self.smis:
				self.all_smis.append(self.smi)

		# draw molecules and save png images for display in structure viewer widget
		# --> there's probably a better way to do this...
		self.draw = 0
		for i, self.smi in enumerate(self.all_smis):
			print(i)
			if i == 0:
				struct_type = 'precursor'
				img_path = 'precursor.png'
				print(struct_type, img_path, self.smi)
			else: struct_type, img_path = 'fragment', 'fragment_%s.png' %self.draw
			print 'bgcolor'
			self.molecule.append(cStructure(struct_type, self.smi, img_path))
			self.template = Chem.MolFromSmiles(self.smi)
			drawOptions = DrawingOptions()
			drawOptions.bgColor = (0,0,0)
			Draw.MolToFile(self.template, img_path, options = drawOptions)

			# append to class level list
			self.structures.append(img_path)

			# append to GUI fragments listbox - file name only
			self.NTPS_frag_list.addItem(QtGui.QListWidgetItem(img_path.split('/')[-1].strip().split('.')[0]))

			self.draw += 1

		self.autocomplete_mod_properties()

		# add image of first structure to graphicsview
		scene = QtGui.QGraphicsScene()
		scene.addPixmap(QtGui.QPixmap(self.structures[0]))
		self.NTPS_structure_view.setScene(scene)
		self.NTPS_structure_view.show()

		return

	def autocomplete_mod_properties(self):
		'''
		Estimate atom ranges and m/z boundaries for NTP search based on input molecule fragments
		'''
		data = element_ranges(self.molecule)

		# add mz window data to mz band field
		mz_band = '%s-%s' %(int(data.min_mz)-1, int(data.max_mz)+1)
		self.NTPS_mz_band_entry.setText(mz_band)

		# add MF details to element table
		self.NTPS_mod_comp_table.setRowCount(len(data.atom_dict))

		i = 0
		for key, value in data.atom_dict.iteritems():
			#print(key, value)
			self.NTPS_mod_comp_table.setItem(i, 0, QtGui.QTableWidgetItem(key))
			#self.NTPS_mod_comp_table.setItem(i, 1, QtGui.QTableWidgetItem(str(value[0])))
			#self.NTPS_mod_comp_table.setItem(i, 2, QtGui.QTableWidgetItem(str(value[1])))

			# center text
			self.NTPS_mod_comp_table.item(i,0).setTextAlignment(QtCore.Qt.AlignCenter)
			#self.NTPS_mod_comp_table.item(i,1).setTextAlignment(QtCore.Qt.AlignCenter)
			#self.NTPS_mod_comp_table.item(i,2).setTextAlignment(QtCore.Qt.AlignCenter)
			i += 1

		# update table measurements to fit data
		self.NTPS_mod_comp_table.resizeRowsToContents()
		self.NTPS_mod_comp_table.resizeColumnsToContents()

		return

	def update_structure_viewer(self, direction):
		'''
		Update structure shown on button click
		'''
		# NTPS_frag_list
		# CC2C=C(C(=O)CN(C)c1c(F)cccc1F)COC2
		row_index = self.NTPS_frag_list.currentRow()
		#row_data = self.NTPS_frag_list.currentItem().text()

		if direction is not None:
			if direction == 'up': new_row = row_index - 1
			elif direction == 'down': new_row = row_index + 1
			self.NTPS_frag_list.setCurrentRow(new_row)

		row_data = self.NTPS_frag_list.currentItem().text()

		scene = QtGui.QGraphicsScene()
		scene.addPixmap(QtGui.QPixmap(row_data))
		self.NTPS_structure_view.setScene(scene)
		self.NTPS_structure_view.show()

		print(row_index, row_data)

		return

	def get_ion_types(self):
		'''
		Get sequence ion types to match
		'''
		ion_types = []
		if self.NTPS_ion_types_a.isChecked() == True: ion_types.append('a')
		if self.NTPS_ion_types_b.isChecked() == True: ion_types.append('b')
		if self.NTPS_ion_types_c.isChecked() == True: ion_types.append('c')
		if self.NTPS_ion_types_x.isChecked() == True: ion_types.append('x')
		if self.NTPS_ion_types_y.isChecked() == True: ion_types.append('y')
		if self.NTPS_ion_types_z.isChecked() == True: ion_types.append('z')
		return tuple(ion_types)

	def get_reactive_residues(self):
		'''
		Get allowed reactive residues
		'''
		react_res = []
		if self.NTPS_reactive_residues_C.isChecked() == True: react_res.append('C')
		if self.NTPS_reactive_residues_Y.isChecked() == True: react_res.append('Y')
		if self.NTPS_reactive_residues_W.isChecked() == True: react_res.append('W')
		if self.NTPS_reactive_residues_H.isChecked() == True: react_res.append('H')
		if self.NTPS_reactive_residues_M.isChecked() == True: react_res.append('M')
		if self.NTPS_reactive_residues_K.isChecked() == True: react_res.append('K')
		if self.NTPS_reactive_residues_R.isChecked() == True: react_res.append('R')
		if self.NTPS_reactive_residues_Q.isChecked() == True: react_res.append('Q')
		return react_res

	def get_atom_range_dict(self, silent = False):

		rowCount = self.NTPS_mod_comp_table.rowCount()

		atom_dict = {}

		for row in xrange(rowCount):
			try:
				element_symbol = str(self.NTPS_mod_comp_table.item(row, 0).text())
				value_string = str(self.NTPS_mod_comp_table.item(row, 1).text())
			except:
				if not silent: print 'Warning - Error reading element table row %s - skipping...' %row
				continue
			a = value_string.replace(' ','').split(',')

			l = []

			for i in a:
				if '-' not in i: l.append(int(i))
				else:
					k1, k2 = i.split('-')
					l += range( int(k1), int(k2) + 1 )

			atom_dict[element_symbol] = l

		# for k,v in atom_dict.iteritems():
		# 	print k,v
		return atom_dict

	def ntps_updateMzBand(self):
		try:
			atom_dict = self.get_atom_range_dict(silent = True)
		except:
			return

		minMassD, maxMassD = {}, {}
		for k,v in atom_dict.iteritems():
			minMassD[k] = min(v)
			maxMassD[k] = max(v)

		minMass = mass.calculate_mass(composition = minMassD)
		maxMass = mass.calculate_mass(composition = maxMassD)

		self.NTPS_mz_band_entry.setText('%s-%s'%(int(minMass), int(maxMass)))
		return

	def ntps_select_output(self):
		NTPS_output = QtGui.QFileDialog.getSaveFileName(self, 'Select Output File', self.fileOpenDialogPath)
		self.fileOpenDialogPath = os.path.dirname(str(NTPS_output))
		self.NTPS_output_file_field.setText(str(NTPS_output))
		return

	def run_non_targeted_search(self):
		'''
		Pull non-targeted search parameters from gui and execute script

		TODO - NEED A CHARGE STATE RANGE SPECIFICATION PARAMETER SOMEWHERE

		'''

		HT_in = self.ntps_htInputFile
		HT_charge = int((self.NTPS_HT_file_charge.text()))
		mascot_in = self.ntps_mascotInputFile

		# MS2 selection parameters

		# get mzBand vlaues sorted lowest first
		mz_band = sorted([float(_) for _ in re.findall(r"[-+]?\d*\.\d+|\d+", str(self.NTPS_mz_band_entry.text()))])

		# get max_RME value
		max_RME = float(self.NTPS_max_RME.text())

		# get HT MS2 correlation boundaries
		HT_ms2_mz_tol = float(str(self.NTPS_HT_MS2_mz_offset.text()))
		HT_ms2_rt_tol = float(str(self.NTPS_HT_MS2_rt_offset.text()))

		# get min HT score
		min_HT_score = int(str(self.NTPS_min_HT_score.text()))

		# get ppm tolerance
		ppm_tol = float(self.NTPS_ppm_tol.text())

		# get match tolerance
		match_tol = float(self.NTPS_ms2_match_tol.text())

		# get user defined reactive residues
		ion_types = self.get_ion_types()

		# get user defined ion types
		reactive_res = self.get_reactive_residues()

		# get atom ranges dictionary
		atom_dict = self.get_atom_range_dict()

		# get output file
		outFile = str(self.NTPS_output_file_field.text())

		# get mascot time unit conversion parameter
		toMins = self.NTPS_toMinutes.isChecked()

		args = {
				'htInput' : [[str(HT_in), HT_charge]],
				'mascotInput' : [str(mascot_in)],
				'mzBand' : mz_band,
				'maxRME' : max_RME,
				'minHTScore' : min_HT_score,
				'ppmTol' : ppm_tol,
				'matchTol' : match_tol,
				'ionTypes' : ion_types,
				'reactiveRes' : reactive_res,
				'HT_ms2_mz_tol' : HT_ms2_mz_tol,
				'HT_ms2_rt_tol' : HT_ms2_rt_tol,
				'atomDict' : atom_dict,
				'outFile' : outFile,
				'toMins' : toMins,
				'molFrags' : self.molecule,
				'drugSMILES' : None,
				'pickle' : True
				}


		run_job(self.runner, self.runner_thread, self.q, NTPI.params, args)
		self.timer.start(1000)
		return

	def check_response(self):
		'''
		On timer tick - check response from NTPS
			- if process returns done, load results into RV GUI tab
		'''
		while not self.q.empty():
			update = self.q.get()
			if update[0] == 'done':
				time.sleep(0.1)
				self.runner.p.terminate()
				self.timer.stop()
				self.runner_thread.exit()
				self.load_RV_results(update)
				print('\nResponse returned to guiprocs')
			else:
				print update
		return

	'''
	RV tab functions
	'''
	def NTPS_RV_gui_connections(self):
		self.NTPS_RV_hits.itemSelectionChanged.connect(lambda: self.NTPS_RV_update_plots(None))
		self.NTPS_RV_up.clicked.connect(lambda: self.NTPS_RV_update_plots('up'))
		self.NTPS_RV_down.clicked.connect(lambda: self.NTPS_RV_update_plots('down'))
		self.NTPS_RV_formulae.itemSelectionChanged.connect(lambda: self.update_RV_structure_viewer(None))
		self.NTPS_RV_accept.clicked.connect(self.add_hit_to_accepted_list)
		self.NTPS_RV_remove.clicked.connect(self.remove_hit_from_accepted_list)
		self.NTPS_RV_accept.clicked.connect(lambda: self.NTPS_RV_update_plots('down'))
		self.NTPS_RV_reject.clicked.connect(lambda: self.NTPS_RV_update_plots('down'))
		self.NTPS_RV_write_report.clicked.connect(self.ntps_RV_write_accepted_report)
		self.NTPS_RV_load_results.clicked.connect(self.loadLocal)
		self.NTPS_RV_reset.clicked.connect(self.reset)
		return

	def loadLocal(self):
		'''
		load a picled list of result objects printed from NTPS
		'''
		import pickle
		NTPS_results_file = str(QtGui.QFileDialog.getOpenFileName(self, 'Select NTPS Results File', self.fileOpenDialogPath))

		with open(r'%s' %NTPS_results_file, 'rb') as f:
			data = pickle.load(f)

		# load into gui
		self.load_RV_results(data)
		return

	def reset(self):

		# clear plots
		self.NTPS_RV_HM.clear()
		self.NTPS_RV_MS.clear()

		for i in reversed(range(self.NTPS_RV_hits.rowCount())):
			self.NTPS_RV_hits.removeRow(i)

		for i in reversed(range(self.NTPS_RV_formulae.rowCount())):
			self.NTPS_RV_formulae.removeRow(i)

		for i in reversed(range(self.NTPS_RV_accepted_hits.rowCount())):
			self.NTPS_RV_accepted_hits.removeRow(i)

		return

	def load_RV_results(self, resp):
		'''
		Parse NTPS search results and load into RV tab

			Attributes of results class are:

				self.pep_sequence = pep_sequence
				self.pep_mz = pep_mz
				self.pep_mod_str = pep_mod_str
				self.ht_hit_mz = ht_hit_mz
				self.mod_mass = mod_mass
				self.score = score
				self.formulae = formulae
				self.index = index
				self.theor_mz = theor_mz
				self.theor_int = theor_int
				self.expt_mz = expt_mz
				self.expt_int = expt_int

				# note
				------------------------------------------
				self.formulae = list of formulas objects

					Attributes of the formulas object are
						formula, mass, ppm

		# need to add more retention time data to this object during NTPS searching
		'''
		try:
			initial_results = resp[1] # search result object
			heatmap_data = resp[2] # list of HT local maxima objects

			# append heatmap data to class level list
			for x in heatmap_data: self.heatmap.append(x)
		except TypeError, e:
			initial_results = resp

		# post processing of results
		# Steps:
		# 	1) Sort hits based on decreasing MS/MS correlation score
		#	2) Calculate RME for each hit for each assigned formula
		#	3) Group hits with the same m/z value - i.e. the same modification --- TODO

		# Sort hits based on correlation score
		initial_results.sort(key = lambda x: float(x.score), reverse = True)
		self.hits = [x for x in initial_results]

		# write results to file
		outFile = str(self.NTPS_output_file_field.text())
		#self.write_results(self.hits, outFile)

		'''
		TODO:
		-------------------------------------------------------
		Move write results function and RME function
		into NTPS script
		'''


		# load all hits into hit list
		for i, hit in enumerate(self.hits):
			hit.count = i
			self.NTPS_RV_hits.insertRow(i)
			self.NTPS_RV_hits.setItem(i, 0, QtGui.QTableWidgetItem(str(hit.count)))
			self.NTPS_RV_hits.setItem(i, 1, QtGui.QTableWidgetItem('%.4f' %float(hit.mod_mass)))

			# centrer text
			self.NTPS_RV_hits.item(i,0).setTextAlignment(QtCore.Qt.AlignCenter)
			self.NTPS_RV_hits.item(i,1).setTextAlignment(QtCore.Qt.AlignCenter)

		# update table measurements to fit data
		self.NTPS_RV_hits.resizeRowsToContents()
		self.NTPS_RV_hits.resizeColumnsToContents()
		self.NTPS_RV_hits.horizontalHeader().setStretchLastSection(True)

		# set first row as highlighted
		self.NTPS_RV_hits.selectRow(0)

		# plot heatmap ---> self.heatmap stores a list of HT objects - attributes are mz, rt, amp, score
		try:
			self.NTPS_RV_HM.clear()
			self.NTPS_RV_HM.plot(x = [x.mz for x in self.heatmap], y = [x.rt for x in self.heatmap], pen = None, symbol = 'o')
		except:
			pass

		# run update_plots to plot EIC and MS
		self.NTPS_RV_update_plots(None)

		# update table measurements to fit data
		self.NTPS_RV_hits.resizeRowsToContents()
		self.NTPS_RV_hits.resizeColumnsToContents()
		return

	def add_hit_to_accepted_list(self):
		'''
		Add selected hit to accepted list
		-- NB: the outcome of NTPS is the molecular formula of the CRM
			--> accepted hit list should contain the hit number of the mascot-HT peptide pair as well as the highlighted formula
			--> need to ensure that an element of the fomula tablewidget is highlighted...
		'''

		# get row number of higlighted peptide hit:
		highlighted_peptide_row = self.NTPS_RV_hits.selectionModel().selectedRows()[0].row()

		# get the index of this hit
		peptide_item = int(str(self.NTPS_RV_hits.item(highlighted_peptide_row, 0).text()))

		# get row number of highlighted formula
		# rows are zero indexed ---> row index should directly correspond to list index in this case
		highlighted_formla_row = int(str(self.NTPS_RV_formulae.selectionModel().selectedRows()[0].row()))

		# get formula data
			# hit index = peptide_item
			# formula index = highlighted_formula_row

		accepted_formula = self.hits[peptide_item].formulae[highlighted_formla_row]

		# get row count of accepted list
		rowCount = self.NTPS_RV_accepted_hits.rowCount()
		self.NTPS_RV_accepted_hits.insertRow(rowCount)

		self.NTPS_RV_accepted_hits.setItem(rowCount, 0, QtGui.QTableWidgetItem(str(peptide_item)))
		self.NTPS_RV_accepted_hits.setItem(rowCount, 1, QtGui.QTableWidgetItem(str(accepted_formula.formula)))
		return

	def remove_hit_from_accepted_list(self):
		# find highlighted row
		highlighted_row = self.NTPS_RV_accepted_hits.selectionModel().selectedRows()[0].row()
		self.NTPS_RV_accepted_hits.removeRow(highlighted_row)
		return

	def write_results(self, results, outFile):

		of1 = open(outFile, 'wt')

		for i, r in enumerate(results):
			of1.write('Summary results for hit %s\n' %i)
			of1.write('====================================\n')
			of1.write('Assigned Mod Mass: %s\n' %r.mod_mass)
			of1.write('Correlation Score: %s\n' %r.score)
			of1.write('Sequence: %s\n' %r.pep_sequence)
			of1.write('Modification Site: %s (%s)\n' %(r.modified_residue, r.modified_residue_index))
			of1.write('HiTIME Score: %s\n' %r.HT_score)
			of1.write('Mascot Score: %s\n' %r.mascot_score)
			of1.write('Mascot Identity/Homology: %s / %s\n' %(r.identity, r.homology))
			of1.write('HiTIME mz: %s\n' %r.ht_hit_mz)
			of1.write('Hitime rt: %s\n' %r.ht_hit_rt)
			of1.write('Mascot mz: %s\n' %r.pep_mz)
			of1.write('Mascot rt: %s\n' %r.pep_rt)

			of1.write('\n')

			of1.write('Possible Formulae\n')
			of1.write('-------------------\n')
			for j, f in enumerate(r.formulae):
				of1.write('Candidate %s: formula: %s, ppm: %s, RME: %s\n' %(j, f.formula, f.ppm, f.RME))

			of1.write('\n\n')

		of1.close()
		return

	def NTPS_RV_update_plots(self, direction):
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

		# TODO
		# --------------------------------------------------------
		# This function is called twice upon every event signal (click, button press etc)
		# Seems to be that this is documented for signals that have optional arguments. See:
		# http://stackoverflow.com/questions/14311578/event-signal-is-emmitted-twice-every-time


		# get selected row from RV_peptides table
		if direction is not None:
			index = self.NTPS_RV_hits.selectionModel().selectedRows()[0].row()
			if direction == 'up':
				self.NTPS_RV_hits.selectRow(index-1)
			elif direction == 'down':
				self.NTPS_RV_hits.selectRow(index+1)

		highlighted_row = self.NTPS_RV_hits.selectionModel().selectedRows()[0].row()

		# get matching hit
		hit = None
		#print('searching for hit matching row:')
		for x in self.hits:
			#print('highlighted_row: %s, hit index: %s' %(highlighted_row, x.count))
			if str(x.count) == str(highlighted_row):
				hit = x
				break

		if hit == None:
			return
		else:
			# update summary text fields
			self.NTPS_RV_mod_mass.setText(str(hit.mod_mass))
			self.NTPS_RV_correlation_score.setText(str(hit.score))
			self.NTPS_RV_sequence.setText(str(hit.pep_sequence))
			self.NTPS_RV_mascot_score.setText(str(hit.mascot_score))
			self.NTPS_RV_HT_score.setText(str(hit.HT_score))
			self.NTPS_RV_mod_site.setText('%s (%s)' %(str(hit.modified_residue), str(hit.modified_residue_index)))
			self.NTPS_RV_mascot_identity_homology.setText('%s/%s' %(str(hit.identity), str(hit.homology)))
			self.NT_RV_delta_RT.setText(str(hit.delta_rt)) #NT_RV_delta_RT
			self.NTPS_RV_HT_MZ.setText(str(hit.ht_hit_mz))
			self.NTPS_RV_HT_RT.setText(str(hit.ht_hit_rt))
			self.NTPS_RV_Mascot_MZ.setText(str(hit.pep_mz))
			self.NTPS_RV_Mascot_RT.setText(str(hit.pep_rt))

			# clear formulae tableWidget
			for i in reversed(range(self.NTPS_RV_formulae.rowCount())):
				self.NTPS_RV_formulae.removeRow(i)

			# add hit formula matches to list
			for i, f in enumerate(hit.formulae):

				self.NTPS_RV_formulae.insertRow(i)
				self.NTPS_RV_formulae.setItem(i, 0, QtGui.QTableWidgetItem(str(f.formula)))
				self.NTPS_RV_formulae.setItem(i, 1, QtGui.QTableWidgetItem(str(int(f.ppm))))
				self.NTPS_RV_formulae.setItem(i, 2, QtGui.QTableWidgetItem(str(int(f.RME))))

				# centrer text
				self.NTPS_RV_formulae.item(i,0).setTextAlignment(QtCore.Qt.AlignCenter)
				self.NTPS_RV_formulae.item(i,1).setTextAlignment(QtCore.Qt.AlignCenter)
				self.NTPS_RV_formulae.item(i,2).setTextAlignment(QtCore.Qt.AlignCenter)

			# set first row of formulae list highlighted
			self.NTPS_RV_formulae.selectRow(0)

			# load first structure into structure viewer
			self.update_RV_structure_viewer(None)

			# update table measurements to fit data
			self.NTPS_RV_formulae.resizeRowsToContents()
			self.NTPS_RV_formulae.resizeColumnsToContents()
			self.NTPS_RV_formulae.horizontalHeader().setStretchLastSection(True)

			# plot heatmap data
			try:
				self.NTPS_RV_HM.clear()
				self.NTPS_RV_HM.plot(x = np.asarray([x.mz for x in self.heatmap]), y = np.asarray([x.rt for x in self.heatmap]), pen = None, symbol = 'o')
			except:
				pass

			# plot MS
			self.NTPS_RV_MS.clear()

			# draw theoretical data
			self.NTPS_RV_MS.plot(x = np.repeat(hit.theor_mz, 3), y = np.dstack((np.zeros(hit.theor_int.shape[0]), hit.theor_int, np.zeros(hit.theor_int.shape[0]))).flatten(), connect='all', pen = self.qp)

			# draw experimental data
			self.NTPS_RV_MS.plot(x = np.repeat(hit.expt_mz, 3), y = np.dstack((np.zeros(hit.expt_int.shape[0]), hit.expt_int, np.zeros(hit.expt_int.shape[0]))).flatten(), connect='all', pen = 'b')

			#self.NTPS_RV_MS.plot(x = np.repeat(hit.expt_mz, 3), y = np.dstack((np.zeros(hit.expt_int.shape[0]), hit.expt_int, np.zeros(hit.expt_int.shape[0]))).flatten(), connect='all', pen = 'b')

		return

	def update_RV_structure_viewer(self, direction):
		'''
		Update structure shown on button click
			- need to retrieve the selected hit first then get the formula list element highlighted

			Attributes of results (hit) object are:

				self.pep_sequence = pep_sequence
				self.pep_mz = pep_mz
				self.pep_mod_str = pep_mod_str
				self.ht_hit_mz = ht_hit_mz
				self.mod_mass = mod_mass
				self.score = score
				self.formulae = formulae
				self.index = index
				self.theor_mz = theor_mz
				self.theor_int = theor_int
				self.expt_mz = expt_mz
				self.expt_int = expt_int

				# note
				------------------------------------------
				self.formulae = list of formulas objects

					Attributes of the formulas object are
						formula, mass, ppm

		# need to add more retention time data to this object during NTPS searching
		'''

		# get index of highlighted row
		row_index = self.NTPS_RV_formulae.currentRow()

		# get index of highlighted row
		highlighted_hit_row = self.NTPS_RV_hits.selectionModel().selectedRows()[0].row()

		# get matching hit
		hit = None
		for x in self.hits:
			if str(x.count) == str(highlighted_hit_row):
				hit = x
				break
		try:
			highlighted_formula_row= self.NTPS_RV_formulae.selectionModel().selectedRows()[0].row()
		except:
			return

		# get formula from
		formula = hit.formulae[highlighted_formula_row]

		scene = QtGui.QGraphicsScene()
		wd = os.getcwd()
		file_list = sorted(os.listdir(wd))

		scene.addPixmap(QtGui.QPixmap(str(formula.RME_frag)))
		self.NTPS_RV_structure.setScene(scene)
		self.NTPS_RV_structure.show()

	def ntps_RV_write_accepted_report(self):

		outfile = str(QtGui.QFileDialog.getSaveFileName(self, 'Select Output File', self.fileOpenDialogPath))
		of1 = open(outfile, 'wt')
		rc = self.NTPS_RV_accepted_hits.rowCount()

		for i in range(rc):
			# get hit number of associated row
			row_hit_number = str(self.NTPS_RV_accepted_hits.item(i,0).text())
			r = None
			for x in self.hits:
				if str(x.count) == str(row_hit_number):
					r = x
					break

			if not r: continue

			of1.write('Summary results for hit %s\n' %i)
			of1.write('====================================\n')
			of1.write('Assigned Mod Mass: %s\n' %r.mod_mass)
			of1.write('Correlation Score: %s\n' %r.score)
			of1.write('Sequence: %s\n' %r.pep_sequence)
			of1.write('Modification Site: %s (%s)\n' %(r.modified_residue, r.modified_residue_index))
			of1.write('HiTIME Score: %s\n' %r.HT_score)
			of1.write('Mascot Score: %s\n' %r.mascot_score)
			of1.write('Mascot Identity/Homology: %s / %s\n' %(r.identity, r.homology))
			of1.write('HiTIME mz: %s\n' %r.ht_hit_mz)
			of1.write('Hitime rt: %s\n' %r.ht_hit_rt)
			of1.write('Mascot mz: %s\n' %r.pep_mz)
			of1.write('Mascot rt: %s\n' %r.pep_rt)
			of1.write('\n')
			of1.write('Possible Formulae\n')
			of1.write('-------------------\n')
			for j, f in enumerate(r.formulae):
				of1.write('Candidate %s: formula: %s, ppm: %s, RME: %s\n' %(j, f.formula, f.ppm, f.RME))

			of1.write('\n\n')
		of1.close()
		return

	def cleanup(self):
		import glob
		files = glob.glob('*.png')
		for filei in files:
			os.remove(filei)
		return
