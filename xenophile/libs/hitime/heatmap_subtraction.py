import sys, os

import os, re, math, argparse
import numpy as np
from xenophile.common import *

# Argument defaults
default_RT_tolerance = 0.4 # min
default_mz_tolerance = 0.2 # m/z
default_score_cutoff = 0

parser = argparse.ArgumentParser(description = 'Subtraction of HiTIME outputs')

parser.add_argument('--inTreatment',
			help = 'name of treatment heatmap input',
			type = str)
parser.add_argument('--inControl',
			help = 'name of control heatmap input',
			type = str)
parser.add_argument('--outFile',
			help = 'name of output file',
			type = str)
parser.add_argument('--rtTolerance',
			help = 'tolerance for retention time differences',
			default = default_RT_tolerance,
			type = float)
parser.add_argument('--mzTolerance',
			help = 'tolerance for MZ differences',
			default = default_mz_tolerance,
			type = float)
parser.add_argument('--scoreCutoff',
			help = 'Scores below this level will be removed prior to subtraction',
			default = default_score_cutoff,
			type = float)


def do_subtraction(treatment, control, options, scoreWise = True):

	# run point subtraction
	Drt = options.rtTolerance
	Dmz = options.mzTolerance

	print Dmz, Drt

	# make list of rts in treatment data
	rts = sorted(list(set(treatment['rt'])))

	subtractedData = []

	for rt in rts:

		# get points at this rt from treatment sample
		tPoints = treatment[np.where(treatment['rt'] == rt)]

		# get points within Drt of rt
		controlIndices = np.where(
									(control['rt'] > (rt-Drt))
									&
									(control['rt'] < (rt + Drt))
								)
		controlPoints = control[controlIndices]

		if len(controlPoints) != 0:
			# iterate through treatment points and get control points w/in Dmz
			for tp in tPoints:
				cps = controlPoints[np.where(
									(controlPoints['mz'] > (tp['mz'] - Dmz))
									&
									(controlPoints['mz'] < (tp['mz'] + Dmz))
									)]

				if len(cps) == 0:
					subtractedData.append(tp)
					#writeData(options, of1, tp) # no control points w/in Dmz

				else:
					# minimise score difference b/w tp and cps
					# - returns cps index of minimised scorre point
					cp_min_index = np.argmin(
											np.absolute(
												cps['score'] - tp['score']
											)
										)
					# get this point
					cp_min_point = cps[cp_min_index]

					# subtract score and write data
					tp['score'] = tp['score'] - cp_min_point['score']
					if scoreWise and tp['score'] > options.scoreCutoff:
						# treatment point overlaps with control point
						# if scoreWise is true, subtract control point score from treatment point score and plot
						print tp
						subtractedData.append(tp)
						#writeData(options, of1, tp)

					else:
						# if scoreWise == False or tp - cp score < 0 - don't write to ouptut
						pass #
		else:
			# no control points found in nearby scans
			# write treatment points to file
			for tp in tPoints:
				subtractedData.append(tp)
				#writeData(options, of1, tp)

	return np.array(subtractedData, dtype = [('mz', '<f8'), ('rt', '<f8'), ('score', '<f8')])

def writeData(options, of1, points):

	for point in points:
		try:
			of1.write('%s, %s, %s\n' %(point['mz'], point['rt'], point['score']))
		except:
			try:
				of1.write('%s, %s, %s\n' %(point.mz, point.rt, point.score))
			except:
				print 'Warning - failed writing point: %s' %point
				pass
	return

def write_headers(options, of1):
	of1.write('# Subtraction_treatment_HT_file::: %s\n' %options.inTreatment)
	of1.write('# Subtraction_control_HT_file::: %s\n' %options.inControl)
	of1.write('# Subtraction_scoreCutoff::: %s\n' %options.scoreCutoff)
	of1.write('# Subtraction_mzTolerande::: %s\n' %options.mzTolerance)
	of1.write('# Subtraction_rtTolerance::: %s\n' %options.rtTolerance)

	# write column headers
	of1.write('## mz, rt, score\n')
	return

def main(options, guiMode = False):

	try:
		controlHeaders, control = read_hitime_files(
												options.inControl,
												scoreCutoff = options.scoreCutoff,
												returnNp = True
												)

		treatmentHeaders, treatment = read_hitime_files(
												options.inTreatment,
												scoreCutoff = options.scoreCutoff,
												returnNp = True
												)

		# are inputs raw score data or list of validated hits?
		if controlHeaders.get('validated', None) and treatmentHeaders.get('validated', None):
			scoreWise = False # do not consider scores when running subtraction - binary point overlap
		else: scoreWise = True
	except:
		pass

	# open output file
	of1 = open(options.outFile,'wt')
	write_headers(options, of1)

	data = do_subtraction(treatment, control, options, scoreWise)


	writeData(options, of1, data)
	of1.close()

	return data


class GuiOptions(object):
	def __init__(self, args):
		for k,v in args.iteritems():
			setattr(self,k,v)

def gui_init(q, args, writeOutput = True):
	options = GuiOptions(args)
	data = main(options, guiMode = True)
	if q:
		q.put('done')
	return data

if __name__ == '__main__':
	options = parser.parse_args()
	main(options)
