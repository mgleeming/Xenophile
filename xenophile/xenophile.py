import sys, os
from common import *
from gui import launch
from PyQt4 import QtCore, QtGui
from multiprocessing import Process, Queue


class LaunchDialog (QtGui.QDialog, launch.Ui_Dialog):

	def __init__(self, parent = None, darkMode = False, runner = None):

		super(LaunchDialog, self).__init__(parent)
		self.setupUi(self)
		self.setWindowTitle('Xenophile')

		self.HT_search_button.clicked.connect(lambda: launch_HT_search(darkMode, runner))
		self.mzML_file_viewer_button.clicked.connect(lambda: launch_MZML_viewer(darkMode))
		self.NTPS_button.clicked.connect(lambda: launch_non_targeted_search(darkMode, runner))
		self.targeted_search_button.clicked.connect(lambda: launch_targeted_search(darkMode, runner))
		return


def launch_non_targeted_search(darkMode, runner):
	from gui.NT_guiprocs import NT_search
	window = NT_search(darkMode = darkMode, runner = runner)
	window.show()
	window.exec_()
	return

def launch_targeted_search(darkMode, runner):
	from gui.TS_guiprocs import T_search
	window = T_search(darkMode = darkMode, runner = runner)
	window.show()
	window.exec_()
	return

def launch_HT_search(darkMode, runner):
	from gui.HT_guiprocs import HT_search
	window = HT_search(darkMode = darkMode, runner = runner)
	window.show()
	window.exec_()
	return

def launch_MZML_viewer(darkMode):
	from gui.mzMLV_guiprocs import mzML_view
	window = mzML_view(darkMode = darkMode)
	window.show()
	window.exec_()
	return

def main():
	''' Init listener thread for multiprocessing tasks'''
	runner_thread = QtCore.QThread()
	runner = Runner(start_signal=runner_thread.started)

	''' launch splash screen '''
	app = QtGui.QApplication(sys.argv)

	darkMode = False


	try:
		if len(sys.argv) > 1 and sys.argv[1] == 'light':
			pass
		else:
			import qdarkstyle as qds
			app.setStyleSheet(qds.load_stylesheet(pyside = False))
			darkMode = True
	except Exception, e:
		print str(e)
		pass

	gui = LaunchDialog(darkMode = darkMode, runner = (runner, runner_thread, Queue()))
	gui.show()
	sys.exit(app.exec_())

import os
import pymzml
import numpy as np
from rtree import index
import matplotlib.pyplot as plt

def read_hitime_files(inFile, returnNp = False, **kwargs):

	htHeaders = dict()

	if '.mzML' in inFile:

		spectra = readSpectra(inFile)

		Times, Mzs, Ints = None, None, None
		for spectrum in spectra:
			times, mzs, ints, lvl = spectrum
			times = np.repeat(times,mzs.shape[0])

			# data = np.vstack((time, mzs, ints))
			if Times is not None:
				Times = np.hstack((Times, times))
				Mzs = np.hstack((Mzs, mzs))
				Ints = np.hstack((Ints, ints))
			else:
				Times = times
				Mzs = mzs
				Ints = ints

		data = [Times, Mzs, Ints]
		npData = np.core.records.fromarrays(
											data,
											names = 'rt, mz, score',
											formats = '<f8, <f8, <f8'
											)
	else:

		f = open(str(inFile),'r')
		# count header rows
		headerRows = 0
		htHeaders = dict()
		while True:
			line = f.readline().strip()
			if line.startswith('#') and not line.startswith('##') and not line.startswith('#-'):
				headerRows += 1 # parameter - value pair
				k,v = line.replace('#','').replace(' ','') .split(':::')
				htHeaders[k] = v
				print k,v

			elif line.startswith('#-'): headerRows += 1 # blank/info only line

			elif line.startswith('##'): # column header line
				names = line.replace('##','').replace(' ','').split(',')
				usecols = (0,1,2)
			else:
				# read until the first non-header line
				if ',' in line:	delimiter = ','
				f.seek(0)
				break

		if headerRows == 0:
			# raw data file
			# - get raw data type

			while True:
				line = f.readline().strip()
				if line != '': break

			# get delimiter
			if ',' in line: delimiter = ','
			else: delimiter = ' '

			data = [x for x in line.replace(',', ' ') if x != '']

			if len(data) == 3:
				htRawDataType = 'cpp'
				names = ['mz', 'rt', 'score']
				usecols = (0,1,2)
			else:
				htRawDataType = 'py'
				names = ['rt', 'mz', 'amp', 'score']
				usecols = (0,1,3)
		# get data and place in numpy array
		f.seek(0)

		# TODO >> assign data types here
		npData = np.genfromtxt(
								f,
								delimiter = delimiter,
								names = names,
								skip_header = headerRows,
								usecols = usecols,
								comments = '#', # skip rows starting with #
							)
		f.close()

	scoreCutoff = kwargs.get('scoreCutoff', None)

	# apply score cutoff
	if scoreCutoff:
		npData  = npData[np.where(npData['score'] > scoreCutoff)]

	# sort entries in order of score - highest first
	npData = npData[npData['score'].argsort()[::-1]]

	# note:
	# to sort lowest first
	# npData[npData['score'].argsort()]

	retrieveTop = kwargs.get('retrieveTop', None)
	if retrieveTop: npData = npData[:retrieveTop]

	# return numpy array if requested
	if returnNp: return htHeaders, npData

	# dictionary specifying hitime_hit attributes other than rt,mz,score
	# these will be applied uniformly to all instances
	addHtAttributes = kwargs.get('addHtAttributes', None)

	# make list of hitime hit objects
	data = []
	for x in npData:
		hit = HtHit(x['mz'], x['rt'], x['score'])

		if addHtAttributes:
			for k,v in addHtAttributes.iteritems(): setattr(hit, k, v)

		data.append(hit)

	# clean up
	del npData

	return htHeaders, data

def loadResults(inFile, returnNp = True, returnHeaders = False):
	headers, data = read_hitime_files(inFile, returnNp = returnNp)
	if returnHeaders:
		return headers, data
	else:
		return data

def summarise(data):
	print 'Number of points: %s' %data.shape[0]
	print 'Average Score: %s' %np.average(data['score'])
	print 'Score Min, Max, Average: %s, %s, %s' %(
								min(data['score']),
								max(data['score']),
								np.average(data['score'])
								)
	print 'm/z Min, Max: %s, %s' %(min(data['mz']), max(data['mz']))
	print 'rt Min, Max: %s, %s' %(min(data['rt']), max(data['rt']))
	print ''
	print '10th percentile: %s' % np.percentile(data['score'], 10)
	print '25th percentile: %s' % np.percentile(data['score'], 25)
	print '50th percentile: %s' % np.percentile(data['score'], 50)
	print '75th percentile: %s' % np.percentile(data['score'], 75)
	print '95th percentile: %s' % np.percentile(data['score'], 95)
	print '99th percentile: %s' % np.percentile(data['score'], 99)
	return

def getPercentile(data, X):
	return np.percentile(data['score'], X)

def getTopX(data, X):
	return data[0:X]

def removeLow(data, minScore):
	return data[np.where(data['score'] > minScore)]

def showScoreHistogram(data, logScale = True, saveFig = False, figName = 'Figure1', dpi = 300):
	bins = int(max(data['score'])) - int(min(data['score'])) + 1
	plt.xlabel('Score')
	plt.ylabel('Frequency')
	plt.title('Score Histogram')
	if logScale:
		plt.yscale('log')
	plt.hist(data['score'], bins = bins, log = True)
	plt.show()
	return

def showHeatMap(data, saveFig = False, figName = 'Figure1', dpi = 300):
	'''
	Plot 2D heatmap of HT results
	Input:
		1) HiTIME results array
	'''
	rdata = data[::-1]

	fig, ax = plt.subplots()
	ax.grid(True)
	ax.set_xlabel('m/z')
	ax.set_ylabel('Retention Time (min)')
	cm = plt.cm.get_cmap('gnuplot_r')
	ax = plt.scatter(rdata['mz'], rdata['rt'], c = rdata['score'], cmap = cm, lw = 0)
	plt.colorbar(ax)

	if not saveFig:
		plt.show()
	else:
		wd = os.getcwd()
		fig.savefig(os.path.join(wd, '%s.png' %figName), dpi = dpi)
	return

def findPeaks(data, Drt = 30, Dmz = 5):
	'''
	Runs local maxima detection to select HiTIME peaks
	Input:
		1) HiTIME results array

	Optional Input:
		1) m/z witdh for detecting local max. Default = 2
		2) rt width for detecting local max. Default = 0.5

	Output:
		1) Structured numpy array containing local maxima points
		   fields are: 'm/z', 'rt', 'score'

	Note: structured array data types
	dtype([('mz', '<f8'), ('rt', '<f8'), ('score', '<f8')])
	'''

	# check for an array
	if not hasattr(data, '__len__'): return None

	data  = data[data['score'].argsort()[::-1]]

	filteredResults = []
	idx = index.Index()
	count = 0
	if Drt > 0 and Dmz > 0:
		for x in data:
			score = x['score']
			mz = x['mz']
			rt = x['rt']
			coord = (rt-Drt, mz-Dmz, rt+Drt, mz+Dmz)
			if idx.count(coord) == 0:
				idx.insert(count, coord)
				filteredResults.append(x)
				count += 1
	return np.array(filteredResults, dtype = [('mz', '<f8'), ('rt', '<f8'), ('score', '<f8')])

class HtHit(object):
	def __init__(self, mz, rt, score):
		self.mz = mz
		self.rt = rt
		self.score = score

def getEIC(data, mzML, mzDelta, mzTolerance = 0.03, msLevel = 1):
	'''
	Get EIC plots for HT targets from provided mzML file and msLevel
	Inputs:
		1) numpy array of HiTIME hits or list of HT Chromatogram objects
		2) mzML file to extract EIC
		3) mzDelta mass spacing between heavy and light
		4) - optional - eicWidth (+/- m/z) used in producting EICs
		5) - optional - mzLevel to use in extraction

	Output:
		1) list of Hitime_EIC objects
			- attributes are: mz, rt, score, EIC_L, EIC_H, RT
	'''

	# create Chromatogram instance for each element of numpy array if needed
	if type(data) != 'list':
		results = [HtHit(d['mz'], d['rt'], d['score']) for d in data]
	else: results = data

	RTs = []
	for r in results:
		# temporarily add these attributes for plotting
		# remove before function return
		r.light_ll = r.mz - mzDelta / 2
		r.light_hl = r.mz + mzDelta / 2
		r.heavy_ll = r.mz - mzDelta / 2
		r.heavy_hl = r.mz + mzDelta / 2
		r.EIC_int_high = []
		r.EIC_int_light = []

	# read spectra and get EIC data
	spectra = __readSpectra(mzML, msLevel = msLevel)
	for n, spectrum in enumerate(spectra):
		time, mzs, ints, lvl = spectrum
		RTs.append(time)

		# test if spectrum corresponds to HT hit point
		for res in results:
			res.EIC_rt.append(float(time))
			res.EIC_int_light.append(0)
			res.EIC_int_heavy.append(0)

		for res in results:
			# get mzs indices within window for light and heavy EICs
			light = np.where((mzs > res.light_ll) & (mzs < res.light_hl))
			heavy = np.where((mzs > res.heavy_ll) & (mzs < res.heavy_hl))

			res.EIC_int_light.append(np.sum(ints[light]))
			res.EIC_int_heavy.append(np.sum(ints[heavy]))

	# EIC plots are held in a nested dictary linked to the EIC attributes
	# Keys of top level dictionary are mzML file names
	# Values are sub dictionariy
	#	Keys/values of sub dictionary are targetMZ, targetRT, mzDelta, RTs, EIC_high, EIC_low
	def addEOC(r, mzML, mzdelta, RTs, addTopLevel = True):
		# add EIC attribute to an instance of the Chromatogram class
		subDict = {
				'targetMZ' : r.mz,
				'targetRT' : r.rt,
				'mzDelta' : mzDelta,
				'RTs' : RTs,
				'EIC_heavy' : r.EIC_int_heavy,
				'EIC_light' : r.EIC_int_light
				}

		# clean temporary attributes
		del r.light_ll
		del r.light_hl
		del r.heavy_ll
		del r.heavy_hl
		del r.EIC_int_heavy
		del r.EIC_int_light

		if addTopLevel:
			r.EIC = {mzML : subDict}
			return r
		else:
			r.EIC[mzML] = subDict
			return r

	for r in results:
		if not hasattr(r, 'EIC'):
			addEIC(r, mzML, mzDelta, RTs)
		else:
			addEIC(r, mzML, mzDelta, RTs, addTopLevel = False)

	return results

def __readSpectra(mzml_file, msLevel):

	msrun = pymzml.run.Reader(str(mzml_file))

	for spectrum in msrun:

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
			mzs = np.array(spectrum.mz, dtype="float32")
			ints = np.array(spectrum.i, dtype ='uint64')
			yield time, mzs, ints, lvl

		except:
			continue

def summariseEICs(result):
	'''
	Show summary of eic traces attached to Chromatogram instance 'result'
	-- print to termina: mzML file origina
						 mzDelta
						 target mz rt
						 num datapoints
	'''
	return

def writeResults(data, outFile):

	of1 = open(outFile,'r')
	of1.write('# m/z, rt, score\n')
	for r in range(data.shape[0]):
		mz, rt, score = data[r]
		of1.write('%s, %s, %s\n'%(mz, rt, score))
	of1.close()
	return

def postProcessHT(
			htFile = None,       mzmlFile = None,
			mzDelta = None,    eicWidth = 0.03,
			scoreCutoff = 0,   outFile = None,
			mzWidth = 0.1,     rtWidth = 0.3,
			peakList = False,  rtExclusion = 0,
			plotEICs = False,  usePeptideIsotopeScaling = False,
			):

	if not all((htFile, outFile)):
		print 'Error - one or more parameters is invalid'
		return None

	if plotEICs and not all((mzDelta, mzmlFile)):
		print 'Error - mzDelta and mzML file needed to plot EICs'
		return None

	import libs.hitime.HT_search_postprocessing as HSP

	args = {
		'htIn' : htFile,
		'mzmlFile' : mzmlFile,
		'mzDelta' : mzDelta,
		'eicWidth' : eicWidth,
		'scoreCutoff' : scoreCutoff,
		'outFile' : outFile,
		'mzWidth' : mzWidth,
		'rtWidth' : rtWidth,
		'peakList' : peakList,
		'rtExclusion' : rtExclusion,
		'plotEICs' : plotEICs,
		'usePeptideIsotopeScaling' : usePeptideIsotopeScaling
		}

	print 'Stargint postprocessing'
	return HSP.guiRun(None, args)

def subtractBackground(
				treatmentData = None,  controlData = None,
				outputFile = None,     rtTolerance = 0.3,
				mzTolerance = 0.1,     scoreCutoff = 0,
				):
	'''
	Subtract data in control array from treatment array
	Input:
		1) treatment array
		2) control array
	Output:
		1) Structured numpy array containing subtracted data points
		   fields are: 'm/z', 'rt', 'score'

	'''

	if not all((treatmentData, controlData, outputFile)):
		print 'Error - specify treatment and control data arrays'
		return None

	import libs.hitime.heatmap_subtraction as hmSub

	args = {
		'inTreatment' : treatmentData,
		'inControl' : controlData,
		'outFile' : outputFile,
		'rtTolerance' : rtTolerance,
		'mzTolerance' : mzTolerance,
		'scoreCutoff' : scoreCutoff
	}

	return hmSub.gui_init(None, args)

def nonTargetedMetID(
				htInput = None, mascotInput = None,
				mzBand = None,  maxRME = 100,
				minHTScore = 0, ppmTol = 100,
				matchTol = 0.5, ionTypes = ('b','y'),
				HT_ms2_mz_tol = 1, HT_ms2_rt_tol = 1,
				atomDict = None, outFile = None,
				toMins = True, drugSMILES = None,
				reactiveRes = ('C','M','W','Y','K')
				):

	if not all((htInput, mascotInput, mzBand,
				atomDict, outFile, drugSMILES)):
		print 'Error - one or more inputs are invalid'
		return None

	import libs.non_targeted_PTM.non_targeted_PTM_identifier as NTPI
	'''
	Note:
	htInput is a nested list of filenames and charges
	[[str(HT_in), HT_charge]
	'''
	args = {
		'htInput' : htInput,
		'mascotInput' : mascotInput,
		'mzBand' : mzBand,
		'maxRME' : maxRME,
		'minHTScore' : minHTScore,
		'ppmTol' : ppmTol,
		'matchTol' : matchTol,
		'ionTypes' : ionTypes,
		'reactiveRes' : reactiveRes,
		'HT_ms2_mz_tol' : HT_ms2_mz_tol,
		'HT_ms2_rt_tol' : HT_ms2_rt_tol,
		'atomDict' : atomDict,
		'outFile' : outFile,
		'toMins' : toMins,
		'drugSMILES' : drugSMILES,  # smiles string for input molecule
		'pickle': True # pickle results object list and save with outfile name
	}

	return NTPI.params(None, args)

def targetedProtID(
				# search options
				htInput = None, mascotInput = None, modIndices = None,
				correlationRT = None, correlationMZ = None,  plotResults = False,
				mzDelta = None, outputFile = None,

				# raw HT parser options
				scoreCutoff = None, rtWidth = None, mzWidth = None,

				# EIC and spectral plotting options
				getIndirect = False, mzmlMS2File = None, mzmlFile = None,
				eicWidth = None,

				# misc
				convertToMinutes = False
			):

	if not all((htInput, mascotInput, outputFile,
	 			modIndices, correlationRT, correlationMZ, mzDelta)):
		print 'Error - one or more inputs are invalid'
		return None

	if plotResults and not all((mzDelta, mzmlFile, eicWidth, mzmlMS2File )):
		print 'Error - one or more inputs are invalid'
		return None

	import libs.targeted_PTM.targeted_PTM as TPM

	args = {
		'outputfile' : outputFile,
		'htInput' : htInput,
		'mascotInput' : mascotInput,
		'getIndirect' : getIndirect,
		'scoreCutoff' : scoreCutoff,
		'rtGap' : rtWidth,
		'mzGap' : mzWidth,
		'modIndices' : modIndices,
		'mzmlFile' : mzmlFile,
		'mzmlMS2File' : mzmlMS2File,
		'plotResults' : plotResults,
		'mzDelta' : mzDelta,
		'eicWidth' : eicWidth,
		'correlationRT': correlationRT,
		'correlationMZ' : correlationMZ,
                'convertToMins' : convertToMinutes
    }

	return TPM.params(None, args)

def runHTSearch(
				mzMLFile,
				mzDelta = None, intensityRatio = 1,
				mzWidth = 150, rtWidth = 17,
				mzSigma = 1.5, rtSigma = 1.5,
				logFile = False, fileFormat = 'mzml',
				noScore = False, removeLow = False,
				ppm = 4
				):

	if not all((mzMLFile, mzDelta)):
		print 'Error - one or more inputs are invalid'
		return None

	import libs.hitime.hitime_methods as HTM
	search_params = {}
	search_params['inputFile'] = str(mzMLFile)
	search_params['mzDelta'] = float(mzDelta)
	search_params['intensityRatio'] = float(intensityRatio)
	search_params['mzWidth'] = float(mzWidth)
	search_params['rtWidth'] = float(rtWidth)
	search_params['mzSigma'] = float(mzSigma)
	search_params['rtSigma'] = float(rtSigma)
	search_params['outputFile'] = str(mzMLFile.split('.')[0]+'.htout')
	search_params['logFile'] = logFile
	search_params['format'] = fileFormat
	search_params['noScore'] = noScore
	search_params['removeLow'] = removeLow
	search_params['ppm'] = 4
	search_params['outDir'] = False
	search_params['minSample'] = rtSigma * rtWidth / 2.355

	HTM.gui_init(None, [search_params])

	return search_params['outputFile']


if __name__ == '__main__':
	main()
