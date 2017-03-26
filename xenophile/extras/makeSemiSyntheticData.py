import os, re, sys, math

import pickle
import numpy as np
import random, argparse
from pyteomics import mass, mgf
from pyteomics import mass, fasta, parser, mgf
from subprocess import Popen, PIPE, STDOUT

import matplotlib.pyplot as plt

try:
	from xenophile.libs.mascotparser import msparser
except ImportError, e:
    try:
        import msparser
    except ImportError, e:
    	print('Msparser import error')
    	print str(e)
    	sys.exit()

print 'imports successful'

parser = argparse.ArgumentParser(description='Create twin-ion MGF entries and mzML data file from CYS/IAA entries.')

# MGF creation options
parser.add_argument('--printVarMods',
                    action='store_true',
                    help='print variable modifications and exit')
parser.add_argument('--mascotFile',
                    help='input mascot data file',
                    type=str)
parser.add_argument('--targetMod',
                    help = 'index of IAA modification',
                    type = str)
parser.add_argument('--outMGF',
                    help = 'output MGF file name',
                    type = str)
parser.add_argument('--exptMGFtarget',
                    help = 'experimental MGF file to which synthetic MS2 data is added',
                    type = str)
parser.add_argument('--outStats',
                    help = 'output summary stats file name',
                    type = str)
parser.add_argument('--keepSyntheticMGF',
                    help = 'Keep MGF file containing MS2 of synthetic peptides',
                    action = 'store_true')

# mzML creation options
parser.add_argument('--mzmlOut',
                    help = 'Output mzML file name',
                    default = 'output.mzML',
                    type = str)
parser.add_argument('--templatemzml',
                    help = 'MzML file to which synthetic data will be added',
                    type = str)
parser.add_argument('--mzDelta',
                    help = 'Mass difference (Da) between heavy and light CRMs',
                    default = 3.01005,
                    type = float)
parser.add_argument('--modMass',
                    help = 'Neutral CRM mass',
                    default = 149.04713,
                    type = float)
parser.add_argument('--minIntensity',
                    help = 'Intensity threshold for synthetic data',
                    default = 5000,
                    type = float)

# misc
parser.add_argument('--rt',
                    help = 'Retention time range (low, high)',
                    nargs = '+',
                    type = str)
parser.add_argument('--mz',
                    help = 'mz range (low, high)',
                    nargs = '+',
                    type = str)


class mascotHit (object): #total_index, prot_index, prot_acc, prot_matches, varmods
    '''
    Storage of mascot peptide assignment data
    '''
    def __init__(self, mz, z, rt, MC, score, sequence, mod_string, query, rank, total_index, prot_index, prot_acc, prot_matches, identity, homology, pep_score):
        self.mz = mz
        self.z = z
        self.rt = rt
        self.MC = MC
        self.score = score
        self.sequence = sequence
        self.mod_string = mod_string
        self.query = query
        self.rank = rank
        self.total_index = total_index
        self.prot_index = prot_index
        self.prot_acc = prot_acc
        self.prot_matches = prot_matches
        self.identity = identity
        self.homology = homology
        self.pep_score = pep_score
        self.lightSN = 0
        self.heavySN = 0
        self.lightsignal = 0
        self.heavysignal = 0
        self.write = True

    def __repr__(self):
        return (16 * '%s, ' %(self.mz, self.z, self.rt, self.MC, self.score, self.sequence, self.mod_string, self.query, self.rank, self.total_index, self.prot_index, self.prot_acc, self.prot_matches, self.identity, self.homology, self.pep_score))

class varMod (object):
    '''
    Storage of variable modifications used in mascot search
    '''
    def __init__(self, modIndex, modName, modDelta, modNeutralLoss):
        self.modIndex = modIndex
        self.modName = modName
        self.modNeutralLoss = modNeutralLoss
        self.modDelta = modDelta

    def __repr__(self):
        return (4 * '%s, ' %(self.modIndex, self.modName, self.modDelta, self.modNeutralLoss))

class PeptideFragments(object):
    def __init__(self, a = None, b = None, c = None, x = None, y = None, z = None, correlation_score = None):
        self.a = a
        self.b = b
        self.c = c
        self.x = x
        self.y = y
        self.z = z

def getIntensityFromMGF(mgfFile):
    data = mgf.read(mgfFile)
    mgfData = np.zeros((100000,3))
    for i, d in enumerate(data):
        if i % 10000 == 0: print 'Processing mgf entry: ', i
        entry = np.array([ float(d['params']['rtinseconds']) , d['params']['pepmass'][0], d['params']['pepmass'][1]])
        #mgfData = np.vstack((mgfData, entry))
        mgfData[i][0] = float(d['params']['rtinseconds'])
        mgfData[i][1] = d['params']['pepmass'][0]
        mgfData[i][2] = d['params']['pepmass'][1]

    return mgfData

def get_peptide_results (resfile, mgfDataArray, options):
    '''
    Retrieve peptide assignments and PTM specifications from mascot .dat file,

    Return values:
        1) list of mascot_hit objects
        2) list of varMod objects
    '''

    # get file header data
    params = resfile.params()

    # get mgf rt vector
    mgfRTs = mgfDataArray[:,0]

    try:
        fixed_mods = params.getMODS()
    except: pass

    # build list of variable modifications and associated mass offsets:
    var_mods = []

    i = 1
    while params.getVarModsName(i):
        modName = params.getVarModsName(i)
        modDelta = params.getVarModsDelta(i)
        modNeutralLoss = params.getVarModsNeutralLoss(i)
        modIndex = i
        var_mods.append(varMod(modIndex, modName, modDelta, modNeutralLoss))
        i += 1

    if options.printVarMods:
        for i in var_mods:
            print i
        sys.exit()

    (scriptName,
     flags,
     minProbability,
     maxHitsToReport,
     ignoreIonsScoreBelow,
     minPepLenInPepSummary,
     usePeptideSummary,
     flags2) = resfile.get_ms_mascotresults_params(msparser.ms_mascotoptions())
    results = msparser.ms_peptidesummary(resfile, flags, 1, 999999999, '', ignoreIonsScoreBelow, minPepLenInPepSummary, '', flags2)
    #results = msparser.ms_peptidesummary(resfile)
    mascot_hits = []

    if usePeptideSummary:
        pepsum = msparser.ms_peptidesummary(
                resfile,                # results file object
                flags,                    # MSRES_group_proteins
                1,                        # minProbability
                999999999,                # maxHits
                '',                        # unigeneIndexFile
                ignoreIonsScoreBelow,     # ignore hits below
                minPepLenInPepSummary,    # minPepLenINPepSummary
                '',                     # singleHit
                flags2)                    # flags2

    total_index = 0

    for x in xrange(1,10000000):
        prot = pepsum.getHit(x)
        # indes, prot_acc, prot_index, prot_matches, varmods
        if prot is not None:
                #print('results for protein hit %x' %x)
                num_peps = prot.getNumPeptides()
                prot_acc = prot.getAccession()
                prot_index = x

                for i in range(1, num_peps + 1):
                    query = prot.getPeptideQuery(i)
                    p = prot.getPeptideP(i)
                    pep = pepsum.getPeptide(query, p)

                    #intensity = resfile.getObservedIntensity(query)

                    if pep.getAnyMatch(): # not sure what this does ---> returns a boolean if any peptide is assigned to this query

                        query = pep.getQuery() # returns index of query
                        queryData = msparser.ms_inputquery(resfile, query)

                        rank = pep.getRank()
                        charge = pep.getCharge()
                        mz = pep.getObserved()
                        seq = pep.getPeptideStr()
                        seq_len = pep.getPeptideLength()
                        score = pep.getIonsScore()
                        #intensity = pep.getTotalIonsIntensity()
                        mod_string = pep.getVarModsStr()
                        prot_matches = pep.getProteins()

                        rt = queryData.getRetentionTimes()

                        miss = pep.getMissedCleavages()

                        identity = results.getPeptideIdentityThreshold(query, 20)
                        homology = results.getHomologyThreshold(query, 20)
                        pep_score = pep.getIonsScore()

                        # TODO: connect UI threshold setting to this conditional
                        if float(score) < float(identity): continue

                        # get 2+ peptides with 1 cysteine
                        if seq.count('C') != 1 or int(charge) != 2: continue

                        # exclude missed cleavages and terminal peptides
                        if seq.count('R') + seq.count('K') != 1: continue

                        # need to count occurrances of IAA/C - make sure CYS is modified w IAA
                        if mod_string.count(str(options.targetMod)) != 1: continue

                        index = np.argmin(np.absolute(mgfRTs - float(rt)))
                        if np.shape < 1:
                            print 'Warning, no MGF intensity found for entry: mz: %s, rt: %s' %(mz, rt)

                        assert float(rt) - float(mgfDataArray[index][0]) < 0.001
                        intensity = mgfDataArray[index][2]

                        if intensity < options.minIntensity: continue

                        hit = mascotHit(float(mz), charge, float(rt), miss, score, seq, mod_string, query, rank, total_index, prot_index, prot_acc, prot_matches, identity, homology, pep_score)

                        hit.exptFragments = []
                        # intensity = 0
                        num_peaks = queryData.getNumberOfPeaks(1)
                        for j in range (1, 1+ num_peaks):
                            peak = [queryData.getPeakMass(1, j), queryData.getPeakIntensity(1, j)]
                            hit.exptFragments.append( peak )
                            # intensity += peak[1]

                        hit.intensity = intensity
                        hit.sequence_mass = float(mass.fast_mass(seq, charge=2))
                        hit.index = total_index
                        mascot_hits.append(hit)
                        total_index += 1


        else: break

    return mascot_hits, var_mods, total_index

def load_mascot_results (if1):
    return msparser.ms_mascotresfile(str(if1))

def get_fragment_mod_masses(frag_mod_str, var_mods):
    '''
    DESCRIPTION
    --------------------------------
    - Determine total mass of all modifications to add to fragment mass
    - Determine mass of candidate CRM in cases where [site] is modified in the native peptide

    INPUTS
    --------------------------------
    site = zero-indexed residue position where candidate CRM is placed
    frag_mod_str = string of integers specifying the indices of the native mods in a mascot-assigned peptide. Integers correspond to modIndex attributes or objects in var_mods list
    var_mods = list of varMods objects. Attributes are: modIndex, modName, modDelta, modNeutralLoss
    ht_mz = m/z value of HT hit. NB - the actual m/z value of the HT hit has been replaced by the precursor m/z value of the nearest MGF query spectrum. The MGF value appears to give better detection of the monoisotopic mass
    pep_mz = the CALCULATED m/z value of the Mascot assigned peptide sequence. No modifications have been incorporated into this number yet.

    RETURNS
    --------------------------------
    1) a list of length (len(frag_mod_str)) that contains the mass adjustments to each residue in a fragment
    2) a float corresponding to the determined CRM mass
    '''

    # determine mass of all native mascot-assigned modifications:
    assigned_mod_mass = 0
    for i in frag_mod_str:
        if i != '0':
            for mod in var_mods:
                if str(mod.modIndex) == i:
                    assigned_mod_mass += mod.modDelta

    # add sum of native mod masses to pep_mz
    #print 'need to replace /2 with charge here'

    # get masses of mods referenced in mod_string
    frag_mod_mass = []
    for i, resMod in enumerate(frag_mod_str):

        # if value in frag_mod_str == 0 ---> No modofication present at this site in Mascot peptide
        if resMod == '0':
            frag_mod_mass.append(float(0))

        # if valeu != 0 ---> Some Mascot-assigned modification is present - find mass of this mod in var_mods list
        else:
            for mod in var_mods:
                if str(mod.modIndex) == resMod:
                    frag_mod_mass.append(mod.modDelta)

    return frag_mod_mass

def get_fragments (peptide, mod_string, var_mods, types = ('b','y'), maxcharge = 2):
    '''
    Generate theoretical sequence ions for a given peptide sequence string
        - types argumnet to be replaced by user specified match ions
        - maxcharge to be replaced by user specified value

    Return Value:
        - List of peptideFragment objects - object for each potential modification site
        - peptideFragment attributes are:

            self.residue = residue
            self.residue_index = residue_index
            self.CRM_mass = CRM_mass
            self.correlation_score = correlation_score

            self.a = a
            self.b = b
            self.c = c
            self.x = x
            self.y = y
            self.z = z

            NB: a,b,c,x,y,z are nested lists - each sublist has the structure ['FRAGMENT SEQUENCE', m/z of fragment]


    NOTE: len(mod_string) == len(peptide) + 2. The two extra entries define modifications at the N and C termini of the pepitde
        ---> for development purposes, remove these

    '''

    mod_string = mod_string[1:len(mod_string)-1]

    #mod_mass = float(ht_hit.mz) - float(peptide.mz) *2  # need to be calculated in the rolling mod function to account for the presence of native pep mods

    assert len(peptide) == len(mod_string)

    # get the mod string for residues in this fragment
    frag_mod_str = list(mod_string)

    # calculate mass of unmodified peptide
    calc_pep_mz = mass.fast_mass(peptide, charge = 2)

    # create a list of masses to add/subtract from each residue
    frag_mod_mass = get_fragment_mod_masses(frag_mod_str, var_mods)

    pepFrags = PeptideFragments()

    a,b,c,x,y,z = [],[],[],[],[],[]

    # generate fragment ions and apply mods (native or CRM) as needed
    for i in xrange(1,len(peptide)):
        for ion_type in types:
            for charge in xrange(1,maxcharge):
                if ion_type in 'abc':

                    # generate pure sequence of fragment
                    frag = peptide[:i]

                    # get mass of relevant mods
                    mods = frag_mod_mass[:i]

                    # calculate mass of base fragment
                    mz = mass.fast_mass(peptide[:i], ion_type = ion_type, charge = charge)

                    # add total mass of modifications - including CRM
                    mz = mz + sum(mods)/charge

                    b.append([peptide[:i], float(mz), ion_type, charge])

                if ion_type in 'xyz':

                    # generate pure sequence of fragment
                    frag = peptide[i:]

                    # get mass of relevant mods
                    mods = frag_mod_mass[i:]

                    # calculate mass of base fragment
                    mz = mass.fast_mass(peptide[i:], ion_type = ion_type, charge = charge)

                    # add total mass of modifications - including CRM
                    mz = mz + sum(mods)

                    y.append([peptide[i:], float(mz), ion_type, charge])


    if len(a) != 0: pepFrags.a = a
    if len(b) != 0: pepFrags.b = b
    if len(c) != 0: pepFrags.c = c
    if len(x) != 0: pepFrags.x = x
    if len(y) != 0: pepFrags.y = y
    if len(z) != 0: pepFrags.z = z

    return pepFrags

def get_modified_fragments(p):

    # Experimental mgf data = p.exptFragments
    # nested list ---> each sublist = [mz, int]

    # theoretical data
    # nested list ---> each sublist = [fragment_sequence, mz, ion type, charge]

    # get all theroetical fragments
    theorPeptide = p.theorFragments

    theorFrags = [theorPeptide.a, theorPeptide.b, theorPeptide.c, theorPeptide.x, theorPeptide.y, theorPeptide.z]
    theorFrags = [x for x in theorFrags if x is not None]
    theorFrags = [i for l in theorFrags for i in l]

    modifiedFragments = []
    num = 0
    # iterate through experimental fragments
    for eF in p.exptFragments:
        eMass, eInt = float(eF[0]), int(eF[1])
        assigned = False

        minError_mgf_entry = None
        minError = None

        # check for cys-containing ions
        for tF in theorFrags:

            if 'C' not in tF[0]: continue

            tMass = float(tF[1])
            if math.fabs(eMass - tMass) < 0.2:
                # Match found - add to assigned ions
                z = int(tF[3])
                mzD = (149.04713 - 57.02146) / z

                modMz = eMass + mzD
                modInt = eInt

                if minError is None and minError_mgf_entry is None:
                    minError = abs(eMass - tMass)
                    minError_mgf_entry = [modMz, modInt] + eF + tF
                elif abs(eMass - tMass) < minError:
                    minError = abs(eMass - tMass)
                    minError_mgf_entry = [modMz, modInt] + eF + tF

        if minError_mgf_entry is not None:
            modifiedFragments.append(minError_mgf_entry)
            assigned = True
            num += 1

        if not assigned: modifiedFragments.append(eF)

    #for x in modifiedFragments: print x
    return modifiedFragments

def make_mgf_entry(modifiedFragments, p):

    # recore rt, mz, int data for newly added peptides
    new_entries = []

    fragMZ = [x[0] for x in modifiedFragments]
    fragINT = [x[1] for x in modifiedFragments]
    pepMass = p.mz + ((149.04713 - 57.02146) / p.z)
    pepRT = float(p.rt) + random.randint(30, 5*60)
    pepHeader = 'synth_%s_original_%s_%s_%s_%s' %(p.index, p.mz, p.rt, p.sequence, p.mod_string)

    assert int(p.z) == 2

    mgfEntry = {
        'm/z array': fragMZ,
        'intensity array': fragINT,
        #'charge array': charges,
        'params': {
            'TITLE': '%s_%s_%s_%s' % (p.index, pepMass, pepRT, pepHeader),
            'PEPMASS': pepMass,
            'RTINSECONDS': pepRT,
            'CHARGE': '2+'
        }
    }

    p.newmz = pepMass
    p.newrt = pepRT

    plot = False

    if plot:

        expMZ, expInt = [],[]
        for eF in p.exptFragments:
            eMass = float(eF[0])
            eInt = int(eF[1])
            expMZ.append(eMass)
            expInt.append(eInt)

        #if len(expMZ) > 0:
        plt.cla()
        plt.clf()
        ml, sl, bl = plt.stem(fragMZ,fragINT, markerfmt = ' ', label = 'expt')
        ml, sl, bl = plt.stem(expMZ,expInt, markerfmt = ' ', label = 'synth')
        plt.setp(sl, color = 'r')
        plt.savefig('figs/peptide_%s_%s_%s_msms.png' %(p.index, p.rt, p.mz))

    return mgfEntry

def writeMGF(spectra, outputFile):

    headers = {
        'COM': 'OpenMS_search',
        'USERNAME': 'OpenMS',
        'FORMAT': 'Mascot generic',
        'TOLU': 'Da',
        'ITOLU': 'Da',
        'FORMVER': '1.01',
        'DB': 'MSDB',
        'SEARCH': 'MIS',
        'REPORT': 'AUTO',
        'CLE': 'Trypsin',
        'MASS': 'monoisotopic',
        'INSTRUMENT': 'Default',
        'PFA': '1',
        'TOL': '3',
        'ITOL': '0.3',
        'TAXONOMY': 'All entries',
        'CHARGE': '1,2,3'
    }

    mgf.write(spectra=spectra, output=outputFile)

    return

def combineMGFfiles(synthMGF, exptMGF):
    of1 = open(synthMGF, 'wt')

    with open(exptMGF,'r') as if1:
        for l in if1:
            of1.write('%s\n' %l.strip())
    of1.write('\n')
    with open(synthMGF + '.synth','r') as if2:
        for l in if2:
            of1.write('%s\n' %l.strip())
    of1.close()
    return

class Spectrum(object):
    def __init__(self, rt, mzs, ints, original):
        self.rt = rt
        self.mzs = mzs
        self.ints = ints
        self.original = original

def getTemplateRts(filei):
    spectra = MZMLtoSpectrum(filei, rtOnly = True)
    return np.array([s.rt for s in spectra])

def prepareSpectra(peptides, modMass, templateFile = None):

    if templateFile:
        templateRTs = getTemplateRts(templateFile)

    # should probably make this a user variable
    # NAPQIp2 = 149.04713 / 2
    # NAPQI13Cp2 = 155.06781 / 2

    for i, p in enumerate(peptides):
        if i % 100 == 0: print 'Preparing Spectrum: ', i
        spectra = MZMLtoSpectrum(p.mzml)
        filerts = []
        fileints = None

        mzOffset = modMass / p.z

        for spectrum in spectra:
            if len(spectrum.mzs) != 0:
                if len(spectrum.mzs) < 30: continue # seems that some short surver

                #p.synthMZs = spectrum.mzs + NAPQIp2
                p.synthMZs = spectrum.mzs + mzOffset

                p.MZmax = np.max(p.synthMZs)
                p.MZmin = np.min(p.synthMZs)

                # get mz of maximal intensity isotope
                maxIntIndex = np.argmax(spectrum.ints)
                p.maxIntMz = p.synthMZs[maxIntIndex]

                filerts.append(spectrum.rt)

                if fileints == None:
                    fileints = spectrum.ints
                else:
                    fileints = np.vstack( (fileints, spectrum.ints) )

                # add in heavy peak later

        # approximate centre of data and rescale retention times
        # peaks may not all have the same number of points
        filerts = range(len(filerts))
        medianRT = len(filerts) / 2
        filerts = np.array(filerts) - medianRT

        # get nearest RT in template file to target
        index = np.argmin(np.absolute(templateRTs - p.newrt))
        # map synthetic rt round this value
        try:
            filerts = templateRTs[filerts + index]
        except:
            filerts = templateRTs[np.arange(10)]
        # normalise intensities and scale to target intensity
        # ---> these are already set at the correct intensity -> specified in FASTA file
        fileints = fileints / np.max(fileints) * p.intensity
        #print np.max(fileints)

        # add to peptide
        p.synthRTs = filerts
        p.RTmax = np.max(p.synthRTs)
        p.RTmin = np.min(p.synthRTs)
        p.synthINTs = fileints
        p.mSig = np.max(p.synthINTs)

    return peptides

def writeCombinedMzML(peptides, outputFile, templateFile, mzDelta, minIntensity):
    import pyopenms, copy
    print 'Writing combined mzML file'

    # generator for template file spectra
    print 'Opening template file'
    spectra = MZMLtoSpectrum(templateFile)

    # output_file = pyopenms.MzMLFile()
    # output_experiment = pyopenms.MSExperiment()

    # can't use usual writer - stores all data in memory before
    # writing to mzML at once ---> too memory intensive

    # use consumer to write spectra on the fly
    consumer = pyopenms.PlainMSDataWritingConsumer(outputFile)

    def replaceData(targetmz, mzs, ints, MZmin, MZmax, synthMZs, synthINTs, minIntensity, peptide):

        # # 1 remove synth data below threshold
        # maparray = np.where(synthINTs > minIntensity)
        # synthMZs = synthMZs[maparray]
        # synthINTs = synthINTs[maparray]

        # check some data is still present
        if synthMZs.shape[0] < 1:
            return mzs, ints

        noise = np.empty(0)
        maxSignal = 0

        maxSN = 0

        # section data into isotopes
        for i in range(20):

            i = float(i)/2
            isotopeMin = MZmin + i - 0.05
            isotopeMax = MZmin + i + 0.3

            if isotopeMin > MZmax: break

            # get synthetic data chuncks
            indexArray = np.where(  (synthMZs >= isotopeMin )
                                    &
                                    (synthMZs < isotopeMax)
                                )

            isotopemzs = synthMZs[indexArray]
            isotopeints = synthINTs[indexArray]

            if isotopemzs.shape[0] < 1:
                #print 'no isotope peaks for peptide %s: mz: %s, rt: %s i-value: %s' %(peptide.index, peptide.newmz, peptide.newrt, i)
                #print 'isotopeMin: %s, isotopeMax: %s' %(isotopeMin, isotopeMax)
                peptide.write = False
                continue

            # get experimental data within isotope boundaries
            minIsotope = np.min(isotopemzs)
            maxIsotope = np.max(isotopemzs)

            indexArray = np.where( (mzs > minIsotope) & (mzs < maxIsotope) )

            isotopeExptlMzs = mzs[indexArray]
            isotopeExptlInts = ints[indexArray]

            # match exptl data to nearest synth mz
            # and increment corresponding xynth int
            for i in range(isotopeExptlMzs.shape[0]):
                eMz = isotopeExptlMzs[i]
                index = np.argmin(np.absolute(synthMZs - eMz))
                synthINTs[index] += isotopeExptlInts[i]


            overlap = np.where( (mzs > minIsotope)
                                &
                                (mzs < maxIsotope)
                                )

            # get overlapping ints
            noise = ints[overlap]

            signal = np.max(isotopeints)

            if signal > maxSignal:
                maxSignal = signal

            mzs = np.delete(mzs, overlap)
            ints = np.delete(ints, overlap)

            mzs = np.concatenate((mzs, synthMZs))
            ints = np.concatenate((ints, synthINTs))

        sortmap = np.argsort(mzs)
        mzs = mzs[sortmap]
        ints = ints[sortmap]

        return mzs, ints, maxSignal

    print 'writing spectra'
    for i, s in enumerate(spectra):

        if i % 100 == 0: print 'Processing spectrum: %s, RT: %.2f' %(i, s.rt/60)

        mzs = s.mzs
        ints = s.ints

        for p in peptides:
            if s.rt > p.RTmax or s.rt < p.RTmin: continue

            for rti, rt in enumerate(p.synthRTs):
                if rt == s.rt: # be careful here

                    synthMZs = p.synthMZs
                    synthINTs1 = p.synthINTs[rti]
                    synthINTs2 = copy.deepcopy(synthINTs1)

                    # add light peak
                    mzs, ints, signal = replaceData(p.newmz, mzs, ints, p.MZmin, p.MZmax, synthMZs, synthINTs1, minIntensity, p)

                    if signal > p.lightsignal:
                        p.lightsignal = signal

                    # add heavy peak
                    mzs, ints, signal = replaceData(p.newmz + mzDelta, mzs, ints, p.MZmin + mzDelta, p.MZmax + mzDelta, synthMZs + mzDelta, synthINTs2, minIntensity, p)

                    if signal > p.heavysignal:
                        p.heavysignal = signal

        new_spectrum = copy.deepcopy(s.original)
        new_spectrum.set_peaks((mzs, ints))

        # use consumer to write new spectrum
        consumer.consumeSpectrum(new_spectrum)

    #output_file.store(outputFile, output_experiment)
    print 'done'
    return peptides

def MZMLtoSpectrum(filename, rtOnly = False):
    spectra = []
    import pyopenms

    # need to use OnDiscMSExperiment to sequentially read
    # spectra without loading whole file into memory
    experiment = pyopenms.OnDiscMSExperiment()
    pyopenms.IndexedMzMLFileLoader().load(filename, experiment)

    for i in range(experiment.getNrSpectra()):
        spectrum = experiment.getSpectrum(i)
        try:
            time = spectrum.getRT()
        except KeyError, e:
            time_prev = time
            if delta_time > 0:
                time += delta_time
            else:
                time += 1.0

        if rtOnly:
            mzData, intData, spectrum = None, None, None

        else:
            (mzData, intData) = spectrum.get_peaks()
            if mzData.shape[0] == 0:
                mzData = np.empty(0, dtype = 'float32')
                intData = np.empty(0, dtype = 'float32')

        yield Spectrum(time, mzData, intData, spectrum)

def main(options):

    '''
    Produce mgf file
    '''
    print '\nGenerating mgf file'
    print '--------------------\n'

    mascotFile = options.mascotFile
    outputMGF = options.outMGF
    statsFileName = options.outStats

    statsFile = open(statsFileName,'wt')

    print 'getting intensities from experimental mgf file'
    mgfDataArray = getIntensityFromMGF(options.exptMGFtarget)

    print 'reading mascot file'
    resfile = load_mascot_results(mascotFile)

    # get assigned peptides that have exactly 1 IAA/Cys
    peptides, varmods, total_peptides = get_peptide_results (resfile, mgfDataArray, options)

    print 'done reading mascot file'

    mgfEntries = []

    for p in peptides:
        # find theoretical fragment ions for each peptide that contain IAA/Cys
        theorFragments = get_fragments(p.sequence, p.mod_string, varmods, types = ('b','y'), maxcharge = 2)
        p.theorFragments = theorFragments
        theorPeptide = p.theorFragments

        # assign cys-containing ions and generate new fragment lists
        modifiedFragments = get_modified_fragments(p)

        mgfEntries.append( make_mgf_entry(modifiedFragments, p) )

    assert len(mgfEntries) == total_peptides

    # write new mgf file
    writeMGF(mgfEntries, outputMGF + '.synth')

    # combine this new MGF file with experimental MGF file
    combineMGFfiles(outputMGF, options.exptMGFtarget)

    print '\nGenerating mzml file'
    print '--------------------\n'

    wd = os.getcwd()
    paramDirectory = os.path.join(sys.path[0], 'params')

    print 'Parsing peptides'
    outputMzML = options.mzmlOut
    iniFile = os.path.join(paramDirectory, 'dataGen.ini')

    print 'Number of peptides is: %s' %len(peptides)

    fastaFiles = []
    mzMLFiles = []

    # create temp file if needed
    tmpFiles = os.path.join(sys.path[0], 'tempFiles')
    if not os.path.exists(tmpFiles):
        os.makedirs(tmpFiles)

    print 'Geneating MS features'
    for pn, peptide in enumerate(peptides):
        if pn % 100 == 0: print 'Generating feature: ', pn

        pepMass = peptide.newmz
        pepRT = peptide.newrt
        pepInt = peptide.intensity

        fasta = os.path.join(tmpFiles, 'temp_%s.fasta' %pn)
        mzml = os.path.join(tmpFiles, 'temp_%s.mzML' %pn)

        peptide.fasta = fasta
        peptide.mzml = mzml

        # write fasta files
        of1 = open(fasta, 'wt')
        of1.write('>seq %s [# intensity= %s, RT=%s #]\n' % (pn, pepInt, pepRT))
        of1.write('%s' % peptide.sequence)
        of1.close()

        # produce mzML file for fasta fragment
        p = Popen( ['MSSimulator' + ' -ini %s' %iniFile + ' -in %s' %fasta + ' -out %s' % mzml], stdout=PIPE, stderr=STDOUT, shell = True, cwd = paramDirectory)
        stdout, nothing = p.communicate()

    # combine mzML file data and write output
    print 'Preparing synthetic spectra'
    peptides = prepareSpectra(peptides, options.modMass, templateFile = options.templatemzml)

    print 'Writing final mzML file'
    outputmzML = os.path.join(wd, outputMzML)
    peptides = writeCombinedMzML(peptides, os.path.join(wd, outputmzML), options.templatemzml, options.mzDelta, options.minIntensity)

    failed = [p for p in peptides if p.write == False]
    passed = [p for p in peptides if p.write == True]

    # write log file
    statsFile.write('# Mascot input file: %s\n' %(mascotFile))
    statsFile.write('# Template mzML input file: %s\n' %(options.templatemzml))
    statsFile.write('# Experimental mgf file: %s\n' %(options.exptMGFtarget))
    statsFile.write('# Output mgf file: %s\n' %(outputMGF))
    statsFile.write('# Output mzML file: %s\n#\n' %(outputmzML))
    statsFile.write('# number of failed peptides is: %s\n#\n' %len(failed))
    statsFile.write('# number of inserted peptides is: %s\n#\n' %len(passed))
    statsFile.write('# new_rt, new_mz, ori_rt, ori_mz, intensity, sequence, new_rt_min, new_rt, max, new_mz_min, new_mz_max, max_intensity_mz, lightsignalmax, heavysignalmax \n')

    for p in peptides:
        if p.write:
            statsFile.write('%.2f, %s, %.2f, %s, %s, %s, %s, %s, %s, %s, %s, %.2f, %.2f\n' %(p.newrt, p.newmz, p.rt, p.mz, p.intensity, p.sequence, p.RTmin, p.RTmax, p.MZmin, p.MZmax, p.maxIntMz, p.lightsignal, p.heavysignal))

    statsFile.write('#\n#\n')

    for p in peptides:
        if not p.write:
            statsFile.write('# %.2f, %s, %.2f, %s, %s, %s, %s, %s, %s, %s, %s, %.2f, %.2f\n' %(p.newrt, p.newmz, p.rt, p.mz, p.intensity, p.sequence, p.RTmin, p.RTmax, p.MZmin, p.MZmax, p.maxIntMz, p.lightsignal, p.heavysignal))

    statsFile.close()

    # clean up
    rubbish = [p.fasta for p in peptides] + [p.mzml for p in peptides]
    if not options.keepSyntheticMGF:
        rubbish.append(outputMGF + '.synth')

    for f in rubbish:
        os.remove(f)
    os.rmdir(tmpFiles)

    return

if __name__ == '__main__':
    options = parser.parse_args()
    sys.exit(main(options))
