import sys

import os, re
from xenophile.common import *

class hit (object):
	def __init__(self):
		return

def getVal(line):
	return re.findall(r"[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", line)

def getString(line):
	return line.split('=')[1].strip()

proteins = []

def parse_file(if1):

	counter = 1
	read_ht_results = False
	results, HT = [], []
	with open(if1,'r') as infile:
		for line in infile:
			if read_ht_results == True and 'END HITIME RESULTS' in line: read_ht_results = False
			if read_ht_results == True:
				rt, mz, score = getVal(line)[:4]
				x = hitime_hit(float(mz), float(rt), float(score))
				HT.append(x)

			if 'Results for Hit' in line:
				results.append(hit())
				results[-1].index = counter
				counter += 1
			elif 'TYPE' in line: results[-1].Type = getString(line)
			elif 'PEPTIDE_MZ' in line: results[-1].mz = getVal(line)[0]
			elif 'PEPTIDE_RT' in line: results[-1].rt = getVal(line)[0]
			elif 'PEPTIDE_CHARGE' in line: results[-1].charge = getVal(line)[0]
			elif 'PEPTIDE_MISSED_CLEAVAGE' in line: results[-1].MC = getVal(line)[0]
			elif 'PEPTIDE_SEQUENCE' in line: results[-1].sequence = getString(line)
			elif 'PEPTIDE_EIC_RT' in line: results[-1].EIC_RT = getVal(line)
			elif 'PEPTIDE_EIC_LIGHT' in line: results[-1].EIC_light = getVal(line)
			elif 'PEPTIDE_EIC_HEAVY' in line: results[-1].EIC_heavy = getVal(line)
			elif 'PEPTIDE_MS_mz' in line: results[-1].pep_MS_mz = getVal(line)
			elif 'PEPTIDE_MS_int' in line: results[-1].pep_MS_int = getVal(line)
			elif 'PEPTIDE_MSN_mz' in line: results[-1].pep_MS2_mz = getVal(line)
			elif 'PEPTIDE_MSN_int' in line: results[-1].pep_MS2_int = getVal(line)
			elif 'HT_MS_mz' in line: results[-1].HT_MS_mz = getVal(line)
			elif 'HT_MS_int' in line: results[-1].HT_MS_int = getVal(line)
			elif 'PROTEIN_ACC' in line: results[-1].prot_acc = getString(line)
			elif 'PROTEIN_INDEX' in line: results[-1].prot_index = getVal(line)[0]
			elif 'PROTEIN_MATCHES' in line: results[-1].prot_matches = getString(line)
			elif 'PEPTIDE_VARMODS' in line: results[-1].mod_string = getString(line)
			elif 'HITIME_MZ' in line: results[-1].HT_mz = getVal(line)[0]
			elif 'HITIME_RT' in line: results[-1].HT_rt = getVal(line)[0]
			elif 'PEPTIDE_IDENTITY_SCORE' in line: results[-1].pep_identity_score = getVal(line)[0]
			elif 'PEPTIDE_HOMOLOGY_SCORE' in line: results[-1].pep_homology_score = getVal(line)[0]
			elif 'PEPTIDE_SCORE' in line: results[-1].pep_score = getVal(line)[0]
			elif 'END RESULTS' in line: read_ht_results = False
			elif 'BEGIN HITIME RESULTS' in line: read_ht_results = True
			else: pass
	return results, HT

def run(ifile = None):
	# testing file

	if ifile is not None:
		results, HT = parse_file(ifile)

	return results, HT
