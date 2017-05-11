import os
import sys
import math
import argparse
import numpy as np
from lxml import etree
from rtree import index
from subprocess import Popen, PIPE
from xenophile.common import *
from PyQt4 import QtCore, QtGui
try:
	from xenophile.libs.mascotparser import msparser
except ImportError, e:
	print('Msparser import error')
	print str(e)


#default full width at half max height of signal for retention time. In number of scans.
DEFAULT_rt_gap = 100  # sec

# default mz range for finding local maxima
DEFAULT_mz_gap = 6 #200

DEFAULT_score_cutoff = 6
DEFAULT_intensity_cutoff = 0

# command line argument parser
parser = argparse.ArgumentParser(description = 'Correlate HiTIME and Mascot outputs to find overlapping peaks')

# File I/O arguments
parser.add_argument('--outputfile',
					type = str,
					help = 'Filename to write output data to')
parser.add_argument('--htInput',
					type = str,
					help = 'HiTIME output file')
parser.add_argument('--mascotInput',
					type = str,
					help = 'Peptide excel output file from mascot')
parser.add_argument('--cleanup',
					action = 'store_true',
					help = 'Remove working files when finished? True/False')
parser.add_argument('--plotData',
					action = 'store_true',
					help = 'Plot full graphical reslts file for each hit')
parser.add_argument('--getIndirect',
					action = 'store_true',
					help = 'Search for indirectly detected modified peptides')

# Hitime output file processing arguments
parser.add_argument('--scoreCutoff',
					type = int,
					action = 'store',
					default = DEFAULT_score_cutoff,
					help = 'Weighted intensity threshold for score filtering')
parser.add_argument('--intensityCutoff',
					type = int,
					action = 'store',
					default = DEFAULT_intensity_cutoff,
					help = 'Absolute intensity threshold for score filtering')
parser.add_argument('--rtGap',
					type = float,
					action = 'store',
					default = DEFAULT_rt_gap,
					help = 'Hitime peak-picking algorithm, rtGap = minimum RT distance between peaks (seconds)')
parser.add_argument('--mzGap',
					type = float,
					action = 'store',
					default = DEFAULT_mz_gap,
					help = 'Hitime peak-picking algorithm, mzGap = minimum mz distance between peaks (mz)')

# EIC arguments - TODO



def infer_HT_hits(HT_data, target_charge):

	inf_data = []

	for datum in HT_data:
		charge_defecit = target_charge - datum.z

		m = datum.mz * datum.z + charge_defecit * 1.00728
		mz = m / target_charge

		hth = hitime_hit(mz, datum.rt, datum.score)
		hth.z = target_charge
		hth.Type = 'light'
		inf_data.append( hth )
	return inf_data

def correlate(options, peptides, HT_data, ofx):
	print('Correlating HiTIME and Mascot results')

	debug = False

	# get correlation data
	Drt = options.correlationRT/2.0
	Dmz = options.correlationMZ/2.0

	NMD = options.mzDelta
	width = options.eicWidth

	results = []

	# this should probably be a user variable
	options.charge_state_range = [2, 3, 4]

	# get charge states represented in HT data
	HT_charges = set([int(x[1]) for x in options.htInput])

	# iterate through charge states and get correlated hits
	for charge in options.charge_state_range:

		print 'Correlating hits for chage state +%s' %charge

		if int(charge) not in HT_charges:
			print 'Inferring HT hit data for +%s' %charge
			inferred_HT_data = infer_HT_hits(HT_data, charge)
		else:
			print 'HT inference not needed - pass'
			inferred_HT_data = None

		HT_charge_data = inferred_HT_data if inferred_HT_data is not None else HT_data

		if debug:
			of1 = open('debug_ht_%s.dat' %charge,'wt')
			of2 = open('debug_mascot_%s.dat' %charge,'wt')

		# init rtree instance
		idxl = index.Index() # 'light' index
		idxh = index.Index() # 'heavy' index

		# HT_data sublists containing data of only the charge of interest
		# ---> add hit points for 'heavy' isotope
		HT_data_charge_light = []
		HT_data_charge_heavy = []

		# print 'HT_data legnth: %s' %len(HT_data_charge_light)
		# print HT_charge_data[0]
		# print type(HT_data[0].z),type(HT_data[0].Type)

		for datum in HT_charge_data:
			if int(datum.z) == int(charge):
				#	print datum
					# append light region
					HT_data_charge_light.append(datum)

					# create heavy region
					heavy_mz = datum.mz + (NMD / datum.z)

					heavy_ht_hit = hitime_hit(heavy_mz, datum.rt,datum.score)
					heavy_ht_hit.z = int(datum.z)
					heavy_ht_hit.Type = 'heavy'
					HT_data_charge_heavy.append(heavy_ht_hit)

					if debug: of1.write('%s, %s\n' %(datum.mz, datum.rt))

		print 'HT_data_charge_light length: %s' %len(HT_data_charge_light)
		print 'HT_data_charge_heavy length: %s' %len(HT_data_charge_heavy)

		# build index of light regions
		for count, datum in enumerate(HT_data_charge_light):
			hrt = datum.rt
			hmz = datum.mz
			coord = (hrt-Drt, hmz-Dmz, hrt+Drt, hmz+Dmz)
			idxl.insert(count, coord)

		# build index of heavy regions
		for count, datum in enumerate(HT_data_charge_heavy):
			hrt = datum.rt
			hmz = datum.mz
			coord = (hrt-Drt, hmz-Dmz, hrt+Drt, hmz+Dmz)
			idxh.insert(count, coord)

		'''
		Find mascot peptides that correlate with HT hits
		'''
		count_lines = 1
		# iterate through mascot peptides
		for peptide in peptides:

			# get ptpdies that
			pmz = float(peptide.mz)
			prt = float(peptide.rt)   #/60
			pcharge = int(peptide.z)

			if pcharge == charge:

				if debug: of2.write('%s, %s\n' %(pmz, prt))

				#print pmz, prt, pcharge
				correlate = (prt - Drt, pmz - Dmz, prt + Drt, pmz + Dmz)

				# returns a list of the indices of the elements that intersect the test box

				'''
				Include function here to check if CRM ID'd is the same (heavy/light) as
				intersection partner
				'''

				hits_light = list(idxl.intersection(correlate))
				hits_heavy = list(idxh.intersection(correlate))

				hits = hits_light + hits_heavy

				if len(hits) > 0:
					#of1.write('0%s %s %s %s\n' %(count_lines, peptide.mz, peptide.rt, peptide.sequence))
					count_lines += 1

					delta = NMD/peptide.z

					# check if any target mods occur in the mod string
					if len(set(str(peptide.mod_string)).intersection(set(str(options.modIndices)))) == 0:
						#print 'no target mods found - moving on'
						continue

					#if peptide.pep_score < peptide.identity:
					#	continue

					# HT_data element that intersected the test peptide
					HT_datum = HT_data_charge_light[hits[0]]

					print '\nPeptide correlation:'
					print pmz, prt, charge, HT_datum.mz, HT_datum.rt, NMD, Drt, Dmz
					#	mz          rt         z ht_mz   ht_rt         NMD    Drt  Dmz
					#	1001.953003 2309.52227 2 752.334 2782.51876851 6.0201 15.0 3.0
					if len(hits) > 1:
						print 'Warning, > 1 ht_data element intersected test peptide'

					res = result('direct',
									peptide.mz,
									peptide.rt,
									HT_datum.mz,
									HT_datum.rt,
									charge,
									peptide.sequence,
									peptide.MC,
									None,
									None,
									None,
									peptide.total_index,
									peptide.prot_index,
									peptide.prot_acc,
									peptide.prot_matches,
									peptide.mod_string,
									peptide.identity,
									peptide.homology,
									peptide.pep_score
									)

					results.append(res)
		if debug:
			of1.close()
			of2.close()
	return results

def mascot_salvage (options, peptides, HT_data, results, ofx):
	# find HT and mascot hist that differ by NAPQI mass
	of1 = open('salvage1.txt','wt')

	# Results Holding List
	for HT_hit in HT_data:
		mz_shift = 57.02146/2 - 149.04713/2
		pep_IAA = HT_hit.mz + mz_shift
		for mascot_hit in peptides:
			if mascot_hit.z == 2:
				if math.fabs(float(pep_IAA) - float(mascot_hit.mz)) < 2:
					if math.fabs(float(HT_hit.rt) - float(mascot_hit.rt)) < 120:

	#					TODO - fix neutral mod delta specification
						NMD = 6.0201
						delta = NMD/mascot_hit.z
						width = 0.03
						ranges = [HT_hit.mz - width, HT_hit.mz + width, HT_hit.mz + delta - width , HT_hit.mz + delta + width]
						res = result('indirect', mascot_hit.mz,mascot_hit.rt, HT_hit.mz, HT_hit.rt, mascot_hit.z, mascot_hit.sequence, mascot_hit.MC, None, ranges, None, None, mascot_hit.total_index, mascot_hit.prot_index, mascot_hit.prot_acc, mascot_hit.prot_matches, mascot_hit.mod_string, mascot_hit.identity, mascot_hit.homology, mascot_hit.pep_score)
						results.append(res)

		#				print('		-- Recovered peak at: %s %s %s %s' %(mascot_hit.mz, mascot_hit.rt, HT_hit.mz, HT_hit.rt))
						of1.write('%s %s %s %s\n' %(mascot_hit.mz ,mascot_hit.rt, HT_hit.mz, HT_hit.rt))

	of1.close()
	return results

# WTF
class result (object):
	def __init__(self,
		res_type,
		pep_mz,
		pep_rt,
		HT_mz,
		HT_rt,
		charge,
		seq,
		MC,
		EIC,
		MS,
		MSn,
		total_index,
		prot_index,
		prot_acc,
		prot_matches,
		mod_string,
		identity,
		homology,
		pep_score):
		self.res_type = res_type
		self.pep_mz = pep_mz
		self.pep_rt = pep_rt
		self.HT_mz = HT_mz
		self.HT_rt = HT_rt
		self.charge = charge
		self.seq = seq
		self.MC = MC
		self.EIC = EIC
		self.MS = MS
		self.MSn = MSn
		self.total_index = total_index
		self.prot_index = prot_index
		self.prot_acc = prot_acc
		self.prot_matches = prot_matches
		self.mod_string = mod_string
		self.identity = identity
		self.homology = homology
		self.pep_score = pep_score

	def __repr__(self):
		return (14 * '%s, ' %(self.res_type, self.pep_mz, self.pep_rt, self.HT_mz, self.HT_rt, self.charge, self.seq, self.MC, self.EIC, self.MS, self.MSn, self.total_index, self.prot_index, self.prot_acc, self.prot_matches, self.mod_string, self.identity, self.homology, self.pep_score))

def get_MS2_and_nearest_MS1(results, options, ms2RTIndex):

	# get indicies of peptide MS2 nad nearest MS1 spectra
	results = resolve_rt_targets_to_spectra(
											results,
											rtIndexArray = ms2RTIndex,
											msLevel = 2,
											mstype = 'MS2'
											)

	# set msLevel to None to bypass msLevel filter
	spectra = readSpectra(options.mzmlMS2File, None)


	print 'running peptide MS2 correlation'
	print 'indices of MS2 target spectra are'
	for res in results:
		print res.pep_rt_index

	count = 0
	# match rt indicies to spectra and save to attribute

	for n, spectrum in enumerate(spectra):
		time, mzs, ints, lvl = spectrum
		count += 1
		for res in results:
			if res.nearest_ms1 == n:
				res.PEP_MS_mz = mzs
				res.PEP_MS_int = ints

			elif res.pep_rt_index == n:

				res.PEP_MS2_mz = mzs
				res.PEP_MS2_int = ints

	#outFile.close()
	print 'number of spectra searched is: %s' %(count)

	return results

class msData(object):
	def __init__(self, rt, mz, intensity, mslevel):
		self.rt = rt
		self.mz = mz
		self.intensity = intensity
		self.mslevel = mslevel

	def __repr__(self):
		return (4 * '%s, ' %(self.rt, self.mz, self.intensity, self.mslevel))

def write_results (results, HT_data, ofx, options):
	print('Writing results file...')
	counter = 1

	print 'length of results list is: %s' %len(results)

	print('		-- Writing peptide results')
	# TODO - write search parameters section
	MS = 1

	# TODO also need to write mods list to results file.
	for result in results:
		if result.res_type == 'direct':
			print >> ofx, 'Results for Hit {}'.format(counter)
			print >> ofx, '..................'.format()
			print >> ofx, 'INDEX = {}'.format(result.total_index)
			print >> ofx, 'TYPE = {}'.format(result.res_type)
			print >> ofx, 'PROTEIN_ACC = {}'.format(result.prot_acc)
			print >> ofx, 'PROTEIN_INDEX = {}'.format(result.prot_index)
			print >> ofx, 'PROTEIN_MATCHES = {}'.format(result.prot_matches)
			print >> ofx, 'PEPTIDE_MZ = {}'.format(result.pep_mz)
			print >> ofx, 'PEPTIDE_RT = {}'.format(result.pep_rt)
			print >> ofx, 'HITIME_MZ = {}'.format(result.HT_mz)
			print >> ofx, 'HITIME_RT = {}'.format(result.HT_rt)
			print >> ofx, 'PEPTIDE_CHARGE = {}'.format(result.charge)
			print >> ofx, 'PEPTIDE_MISSED_CLEAVAGE = {}'.format(result.MC)
			print >> ofx, 'PEPTIDE_SEQUENCE = {}'.format(result.seq)
			print >> ofx, 'PEPTIDE_VARMODS = {}'.format(result.mod_string)
			print >> ofx, 'PEPTIDE_IDENTITY_SCORE = {}'.format(result.identity)
			print >> ofx, 'PEPTIDE_HOMOLOGY_SCORE = {}'.format(result.homology)
			print >> ofx, 'PEPTIDE_SCORE = {}'.format(result.pep_score)

			if options.plotResults:
				print >> ofx, 'PEPTIDE_EIC_RT = {}'.format(result.EIC_rt)
				print >> ofx, 'PEPTIDE_EIC_LIGHT = {}'.format(result.EIC_int_light)
				print >> ofx, 'PEPTIDE_EIC_HEAVY = {}'.format(result.EIC_int_heavy)
				print >> ofx, 'PEPTIDE_MS_mz = {}'.format(list(result.PEP_MS_mz))
				print >> ofx, 'PEPTIDE_MS_int = {}'.format(list(result.PEP_MS_int))
				print >> ofx, 'PEPTIDE_MSN_mz = {}'.format(list(result.PEP_MS2_mz))
				print >> ofx, 'PEPTIDE_MSN_int = {}'.format(list(result.PEP_MS2_int))
				print >> ofx, 'HT_MS_mz = {}'.format(list(result.HT_MS_mz))
				print >> ofx, 'HT_MS_int = {}'.format(list(result.HT_MS_int))

			print >> ofx, 'END RESULTS'.format()
			print >> ofx, ''.format()
			counter += 1

		# elif result.res_type == 'indirect':
		# 	print >> ofx, 'Results for Hit {}'.format(counter)
		# 	print >> ofx, '..................'.format()
		# 	print >> ofx, 'INDEX = {}'.format(result.total_index)
		# 	print >> ofx, 'TYPE = {}'.format(result.res_type)
		# 	print >> ofx, 'PROTEIN_ACC = {}'.format(result.prot_acc)
		# 	print >> ofx, 'PROTEIN_INDEX = {}'.format(result.prot_index)
		# 	print >> ofx, 'PROTEIN_MATCHES = {}'.format(result.prot_matches)
		# 	print >> ofx, 'PEPTIDE_MZ = {}'.format(result.pep_mz)
		# 	print >> ofx, 'PEPTIDE_RT = {}'.format(result.pep_rt)
		# 	print >> ofx, 'HITIME_MZ = {}'.format(result.HT_mz)
		# 	print >> ofx, 'HITIME_RT = {}'.format(result.HT_rt)
		# 	print >> ofx, 'PEPTIDE_CHARGE = {}'.format(result.charge)
		# 	print >> ofx, 'PEPTIDE_MISSED_CLEAVAGE = {}'.format(result.MC)
		# 	print >> ofx, 'PEPTIDE_SEQUENCE = {}'.format(result.seq)
		# 	print >> ofx, 'PEPTIDE_VARMODS = {}'.format(result.mod_string)
		# 	print >> ofx, 'PEPTIDE_EIC = {}'.format(result.EIC)
		# 	print >> ofx, 'PEPTIDE_MS = {}'.format(result.MS)
		# 	print >> ofx, 'PEPTIDE_MSN = {}'.format(result.MSn)
		# 	print >> ofx, 'END RESULTS'.format()
		# 	print >> ofx, ''.format()
		# 	counter += 1

	print('		-- Writing HiTIME summary results')
	print >> ofx, 'BEGIN HITIME RESULTS'
	for hit in HT_data:
		print >> ofx, '{}, {}, {}'.format(hit.rt, hit.mz, hit.score)
	print >> ofx, 'END HITIME RESULTS'
	print('Done!')
	return

def main(options, guiMode = False, queue = None):

	if queue:
		msg = 'Begin targeted search'
		queue.put(msg)

	# init output file
	ofx = open(options.outputfile,'w')

	# Extract HT results
	print 'Extracting HT data'

	HT_data = []

	# iteratively parse HT files
	for htInput in options.htInput:

		htInputFile = htInput[0] # target file name
		htInputCharge = htInput[1] # target file charge

		additional_HT_attributes = dict()
		additional_HT_attributes['z'] = int(htInputCharge)
		additional_HT_attributes['Type'] = 'light'

		headers, this_HT_data = read_hitime_files(
											htInputFile,
											addHtAttributes = additional_HT_attributes
											)

		# TODO: rtWidth, mzWidth and scoreCutoff are only applied in cases
		# where raw HT data is supplied. i.e. Where validated HT hit lists
		# are provided, these parameters are redundant.
		# -- This is not obvious from the GUI.
		# Should parse the HT file headers in the TS_guiprocs and grey-out
		# mzWidth, rtWidth and scoreCutoff parameters in cases where validated
		# peak lists are given.


		# check if this ht data is a verified peak list - default to false if not
		peakList = headers.get('validated', False)

		if not peakList:
			# run automated peak picking if not a peak list
			this_HT_data = get_HT_regions(this_HT_data,
									options.rtWidth, # Drt --- Not sure if these
									0.1, # Dmz --- should be user variables
									minScore = options.scoreCutoff
									)

		HT_data += this_HT_data

	# load mascot file and build RT index
	print 'Extracting Mascot data'
	resfile = load_mascot_results(options.mascotInput)
	#rt_data = get_rt_data(resfile, convertToMinutes = True)

	# get hitime input files
	peptides, var_mods = get_peptide_results(resfile, convertToMinutes = False)

	# run HT/peptide correlation
	print 'Running correlation'
	results = correlate(options, peptides, HT_data, ofx)

	print '%s hits identified by direct correlation' %len(results)
	#results = mascot_salvage(options, peptides, HT_data, results, ofx)


	if options.plotResults and options.mzmlFile:

		ms1RTIndex, unit = build_rt_index(options.mzmlFile)

		# ofx = open('ms1RTIndex.dat','wt')

		# for i in ms1RTIndex:
		# 	ofx.write('%s, %s\n' %(i[0],i[1]))
		# ofx.close()

		#if testing is not True:
		if options.mzmlMS2File:
			ms2RTIndex, unit = build_rt_index(options.mzmlMS2File)

		# ofx = open('ms2RTIndex.dat','wt')
		# for i in ms2RTIndex:
		# 	ofx.write('%s, %s\n' %(i[0],i[1]))
		# ofx.close()

		print 'ms1RTIndex size is: %s' %ms1RTIndex.size
		print 'ms2RTIndex size is: %s' %ms2RTIndex.size

		print 'Resolving MS1 targets'
		results = resolve_rt_targets_to_spectra(
												results,
												rtIndexArray = ms1RTIndex,
												msLevel = 1,
												mstype = 'MS1'
												)

		'''
			1) compare HT hits to rts in ms1 rt index
			2) for ms1 spectra, find ms closes
		'''

		# get EIC and MS1 spectra for HT hits
		results = plot_EIC(results, options)


		print 'Resolving MS2 spectra'
		# get MS2 spectra and nearest MS1 spectra for peptide hits
		results = get_MS2_and_nearest_MS1(results, options, ms2RTIndex)

	print 'Writing results to output file'
	write_results(results, HT_data, ofx, options)
	ofx.close()

	if queue:
		msg = 'Finished targeted search'
		queue.put(msg)
		time.sleep(1)
		queue.put('done')
	return

class GuiOptions(object):
	''' manually create options object when script called from GUI '''
	def __init__(self, args):
		for k,v in args.iteritems():
		#	print k,v
			setattr(self, k, v)

def params(q, args):
	'''
	Init function called from GUI

	data = dictionary of input arguments/values
	'''

	options = GuiOptions(args)
	main(options, guiMode = True, queue = q)

if __name__ == '__main__':
	options = parser.parse_args()
	main(options, guiMode = False)
