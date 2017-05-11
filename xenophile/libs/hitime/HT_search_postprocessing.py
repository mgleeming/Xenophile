import os, re, sys, operator
from datetime import datetime
import numpy as np
import argparse

from xenophile.common import *

parser = argparse.ArgumentParser(description = 'find hits in a HiTIME reuslts file')

parser.add_argument('-htIn',
					help = 'HiTIME results file or peak list for processing',
					required = True,
					type = str)
parser.add_argument('-outFile',
					help = 'Output file name',
					type = str)
parser.add_argument('--mzmlFile',
					help = 'mzML file used in HiTIME searching',
					type = str)
parser.add_argument('--peakList',
					help = 'Input file provided is a list of HT peaks',
					action = 'store_true')
parser.add_argument('--mzDelta',
					help = 'mass difference between heavy and light isotopes. Default = 6.0201',
					default = 6.0201,
					type = float)
parser.add_argument('--rtWidth',
					help = 'retention time width for local maxima detection',
					type = float)
parser.add_argument('--mzWidth',
					help = 'm/z width for local maxima detection',
					type = float)
parser.add_argument('--eicWidth',
					help = 'm/z width used for plotting Extracted Ion Chromatograms. Default = 0.03',
					default = 0.03,
					type = float)
parser.add_argument('--scoreCutoff',
					help = 'HiTIME score threshold. Default = 0',
					default = 0,
					type = float)
parser.add_argument('--rtExclusion',
					help = 'Exclude points for specified time following accepted hit',
					default = 0,
					type = float)


def write_headers(options, unit, htHeaders, of1):
	of1.write('#- BEGIN_POSTPROCESSING_PARAMETERS\n')
	of1.write('# htFile::: %s\n' %options.htIn)
	of1.write('# mzmlFile::: %s\n' %options.mzmlFile)
	of1.write('# mzWidth::: %s\n' %options.mzWidth)
	of1.write('# rtWidth::: %s\n' %options.rtWidth)
	of1.write('# scoreCutoff::: %s\n' %options.scoreCutoff)
	of1.write('# mzDelta::: %s\n' %options.mzDelta)
	of1.write('# eicWidth::: %s\n' %options.eicWidth)
	of1.write('# mzML_time_unit::: %s\n' %unit)
	for k,v in htHeaders.iteritems():
		of1.write('%s::: %s\n' %(k,v))
	of1.write('#- END_POSTPROCESSING_PARAMETERS\n')
	return

def write_ht_hits_to_file(results, of1):
	of1.write("#- BEGIN_RAW_HITIME_RESULTS\n")
	of1.write('## rt, mz, score\n')
	for i in results:
		of1.write('%s, %s, %s\n' %(i.rt, i.mz, i.score))
	of1.write("#- END_RAW_HITIME_RESULTS\n\n")
	return

def write_data_to_file(results, of1):
	of1.write('BEGIN_HITIME_HIT_RESULTS\n')
	for index, i in enumerate(results):
			of1.write('\nBEGIN_Results_for_hit %s\n' %index)
			of1.write('<RT> %.2f\n' %i.rt)
			of1.write('<MZ> %s\n' %i.mz)
			of1.write('<SCORE> %.2f\n' %i.score)
			of1.write('<EIC_RT> %s\n' %['%.3f' %x for x in i.EIC_rt])
			of1.write('<EIC_int_light> %s\n' %i.EIC_int_light)
			of1.write('<EIC_int_heavy> %s\n' %i.EIC_int_heavy)
			of1.write('<MS_mz> %s\n' %[('%.4f' %x) for x in i.HT_MS_mz])
			of1.write('<MS_int> %s\n' %[int(x) for x in i.HT_MS_int])
			of1.write('<END_Results_for_hit> %s\n' %index)
	return

def filter_isotopes(results, options):

	# sort results list in order of increasing m/z value
	results.sort(key = lambda d: d.mz)

	print 'Filtering isotopic hits'
	print 'Length of initial results list: %s' %len(results)

	idx = index.Index()

	Drt = options.rtWidth
	Dmz = options.mzWidth
	rtExclusion = options.rtExclusion
	mzIso = 2.5 # make this a user variable?

	filtered_results = []
	index_counter = 0
	for res in results:
		mz = res.mz
		rt = res.rt

		coord = (rt - Drt, mz - Dmz, rt + Drt, mz + mzIso)

		intersections = list(idx.intersection(coord))

		if len(intersections) == 0:
			idx.insert(index_counter, coord)
			filtered_results.append(res)
			index_counter += 1

		else:
			conflictPoint = filtered_results[intersections[0]]
			conf_mz = conflictPoint.mz
			conf_rt = conflictPoint.rt

			if options.usePeptideIsotopeScaling:
				conflictPoint.score *= res.score

			print 'Isotope exclusion of m/z %s, rt %s - conflict with m/z %s, rt %s' %(mz, rt, conf_mz, conf_rt)

	# re-sort new results list in order of decreasing score
	filtered_results.sort(key = lambda d: d.score, reverse = True)

	print 'Length of filtered results list: %s' %len(filtered_results)

	return filtered_results

def write_accepted_to_file(headers, data, options, of1):

	of1.write('#- Validated HiTIME hit list:\n')

	# write postprocessing headers to results file
	of1.write('#- Postprocessing parameters\n')
	for k, v in headers.iteritems():
		of1.write('# %s::: %s\n'%(k,v))

	of1.write('# htFile::: %s\n' %options.htIn)
	of1.write('# mzmlFile::: %s\n' %options.mzmlFile)
	of1.write('# mzWidth::: %s\n' %options.mzWidth)
	of1.write('# rtWidth::: %s\n' %options.rtWidth)
	of1.write('# scoreCutoff::: %s\n' %options.scoreCutoff)
	of1.write('# mzDelta::: %s\n' %options.mzDelta)
	of1.write('# eicWidth::: %s\n' %options.eicWidth)
	#of1.write('# mzML_time_unit::: %s\n' %unit)

	of1.write('# validated::: 1\n') # indicates that these results have been user-checked

	of1.write('#-\n')

	of1.write('## mz, rt, score\n')
	# get data from rows of accepted hit table
	for hit in data:
		rt = hit.rt
		mz = hit.mz
		score = hit.score
		of1.write('%s, %s, %s\n' %(mz, rt, score))

	return

#@profile
def main(options):

	# Create output file
	of1 = open(options.outFile, 'wt')

	# parse raw HT datafile
	#results = sort_HT_results(options, of1)
	print 'Reading HiITME file'
	htHeaders, results = read_hitime_files(
											options.htIn, # input ht file
											mzWidth = options.mzWidth,
											rtWidth = options.rtWidth,
											scoreCutoff = options.scoreCutoff
											)

	print 'Finished reading HiTIME file'

	# get local maxima for HT search if needed
	if not options.peakList:
		print 'Filtering Isotopes'
		results = get_HT_regions(
								results,
								options.rtWidth,
								options.mzWidth,
								rtExclusion = options.rtExclusion
								)

		# add found local max to headers
		htHeaders['foundLocalMax'] = 1

	print '%s HiTIME hits detected' %len(results)

	# exclude A+1 and A+2 isotopic peaks
	results = filter_isotopes(results, options)

	if options.plotEICs:
		# build index of rt's for all spectra in supplied mzML file
		# unit returns the time value units from the mzML file
		ms1RTIndex, unit = build_rt_index(options.mzmlFile)

		print 'Building RT index'
		results = resolve_rt_targets_to_spectra(
												results,
												rtIndexArray = ms1RTIndex,
												msLevel = 1,
												mstype = 'MS1'
												)

		print 'Plotting EIC'
		# extract EIC and MS spectrumfor each hit
		results = plot_EIC(results, options)

	else:
		unit = None

	print 'Writing data to file'
	# write data to file

	if options.plotEICs:
		# write headers
		write_headers(options, unit, htHeaders, of1)

		# write HT htis to file
		write_ht_hits_to_file(results, of1)

		# write hit MS and EIC data to file
		write_data_to_file(results, of1)

	else:
		write_accepted_to_file(htHeaders, results, options, of1)
	of1.close()

	print 'Finished postprocessing'
	return

#run('rat.max.threshold.1000.m300.scores.sorted.top.100.txt', '../rat.threshold.1000.mzML', 'speed-test.dat')

class GuiOptions(object):
	def __init__(self, args):
		for k,v in args.iteritems():
			setattr(self, k, v)

def guiRun(q, args):
	options = GuiOptions(args)
	main(options)
	if q:
		q.put('done')
	return

if __name__ == '__main__':
	options = parser.parse_args()
	sys.exit(main(options))
