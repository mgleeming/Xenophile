
import numpy as np
class hitime_hit (object):
	def __init__(self, mz, rt, score, amp = None):
		self.mz = mz
		self.rt = rt
		self.amp = amp
		self.score = score
	
	def __repr__(self):
		return ( 4 * '%s, ' %(self.mz, self.rt, self.amp, self.score))

def read_HT_data(inFile, returnNp = False, **kwargs):
	# get hitime file type
	# options are:
	#	- raw data python
	# 	- raw data cpp
	# 	- subtracted data
	# 	- validated hist from gui
	# 	- invalid

	'''
	raw data python
		4 columns	- rt, mz, amp, score

	raw data cpp 
		3 columns	- mz, rt, score

	validated data
		headers 
			- begin with '#'
			- 
		data formatted as mz, rt, score

	subtracted data
		headers - 


	general approach
	----------------------------------

	use numpy np.getfromtxt function
		has arguments: delimiter, names (column names), skip_header .... + more


	if headers encountered in first x lines
		- use this infor to parameterise numpy matrix
		- count number of header lines

	read data into numpy array

	sort by score

	apply score cutoff

	create hitime_hit objects if needed

	'''

	f = open(inFile,'r')


	# count header rows

	headerRows = 0
	htHeaders = dict()
	while True:
		line = f.readline().strip()
		if line.startswith('#'):
			headerRows += 1
			k,v = line.split(':::')
			htHeaders[k] = v

		elif line == '': continue
		else:
			f.seek(0)
			delimiter = ','
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
		else:
			htRawDataType = 'py'
			names = ['rt', 'mz', 'amp', 'score']

		if 'names' not in htHeaders:
			htHeaders['colNames'] = names

	print htHeaders['colNames']

	# get data and place in numpy array
	f.seek(0)
	npData = np.genfromtxt(
							f,
							delimiter = delimiter,
							names = htHeaders['colNames'],
							skip_header = headerRows
						)

	scoreCutoff = kwargs.get('scoreCutoff', None)

	# apply score cutoff
	if scoreCutoff:
		npData  = npData[np.where(npData['score'] > scoreCutoff)]

	# sort entries in order of score - highest first
	npData = npData[npData['score'].argsort()[::-1]]
	
	# note:
	# to sort lowest first 
	# npData[npData['score'].argsort()]

	# return numpy array if requested
	if returnNp: return htHeaders, npData

	# dictionary specifying hitime_hit attributes other than rt,mz,score
	# these will be applied uniformly to all instances
	addHtAttributes = kwargs.get('addHtAttributes', None)

	# make list of hitime hit objects
	data = []
	for x in npData:
		hit = hitime_hit(x['mz'], x['rt'], x['score'])

		if addHtAttributes:
			for k,v in addHtAttributes.iteritems(): setattr(hit, k, v)

		data.append(hit)

	getLocalMax = kwargs.get('getLocalMax', None)	

	if getLocalMax:
		mzWidth = kwargs.get('mzWidth', None)
		rtWidth = kwargs.get('rtWidth', None)
		rtExclusion = kwargs.get('rtExclusion', 0)

		if mzWidth and rtWidth:
			data = getHtRegions(data, rtWidth, mzWidth, rtExclusion)

	# clean up
	del npData
	f.close()

	return htHeaders, data

def getHtRegions(HT_data, Drt, Dmz, rtExclusion = 0):
	if Drt is None or Dmz is None:
		return HT_data

	count = 0
	HT_regions = []

	idx = index.Index()

	for x in HT_data:
		mz = x.mz
		rt = x.rt
		score = x.score
		coord = (rt-Drt, mz-Dmz, rt+Drt+rtExclusion, mz+Dmz)
		if idx.count(coord) == 0:
			idx.insert(count, coord)
			HT_regions.append(x)
			count += 1

	return HT_regions
