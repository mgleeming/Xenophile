import sys
import os, re, math, random
import argparse
import itertools
import resource
import numpy as np
import resource
from rtree import index
import matplotlib.pyplot as plt
from pyteomics import mgf, pepxml, mass
from datetime import datetime

from xenophile.common import *

try:
    from xenophile.libs.mascotparser import msparser
except ImportError, e:
    print('Msparser import error')
    print str(e)
    sys.exit()

# set command line argumentsf
parser = argparse.ArgumentParser(description = 'Non-targeted search for reactive metabolite-peptide adducts')

parser.add_argument('--htInput', type = str, help = 'HiTIME output file')
parser.add_argument('--mascotInput', type = str, help = 'Mascot output .dat file')
parser.add_argument('--mzBand', type = int, default = None, help = 'm/z region relative to a HiTIME hit to search for potential unmodified peptides')
parser.add_argument('--rtBand', type = int, default = None, help = 'Retention time band relative to a HiTIME hit to search for potential unmodified peptides')
parser.add_argument('--minHTScore', type = float, default = None, help = 'Minimum HiTIME score required for a peptide to be considered')
parser.add_argument('--ppmTol', type = float, default = None, help = 'ppm error curoff for matching potential CRM formulae to a mass')
parser.add_argument('--matchTol', type = float, default = None, help = 'm/z tolerance for matching predicted and observed MS2 peaks')
parser.add_argument('--ionTypes', type = str, default = None, nargs = '+', help = 'Ion types to be matched against experimental spectra. Valid options are: a, b, c, x, y, z. Enter each desired ion type separated by a space.')
parser.add_argument('--reactiveRes', type = str, default = None, nargs = '+', help = "Peptide residues that are considereed 'reactive'. Peptides that do not contain at least one of these will not be considered as match candidates. Valid options are: C, Y, W, H, M, K, R, Q. Enter each desired residue separated by a space")
parser.add_argument('--HT_ms2_mz_tol', type = float, default = 0.5, help = 'm/z tolerance for matching HT hits with associated MS2 spectra')
parser.add_argument('--HT_ms2_rt_tol', type = float, default = 100, help = 'rt tolerance for matching HT hits with associated MS2 spectra')

class mgfData (object):
    # class for mgf query data
    def __init__(self, index, title, mzMin, mzMax, intMin, intMax, numVals, numUsed, ms2mz, ms2int, rt, mz, charge):
        self.index = index
        self.title = title
        self.mzMin = mzMin
        self.mzMax = mzMax
        self.intMin = intMin
        self.intMax = intMax
        self.numVals = numVals
        self.numUsed = numUsed
        self.ms2mz = ms2mz
        self.ms2int = ms2int
        self.rt = rt
        self.mz = mz
        self.charge = charge

    def __repr__(self):
        return (13*'%s, ' %(self.index, self.title, self.mzMin, self.mzMax, self.intMin, self.intMax, self.numVals, self.numUsed, self.ms2mz, self.ms2int, self.rt, self.mz, self.charge))

def get_mgf_data(resfiles, convertToMinutes = False):

    # for testing only - load 2 resfiles - one with mods and one without
    query_data = []

    for resfilei in list(resfiles):
        resfile = msparser.ms_mascotresfile(resfilei) # file name = user variable

        for i in xrange(1, 1 + resfile.getNumQueries()):
            q = msparser.ms_inputquery(resfile, i)

            title = q.getStringTitle(i)

            mzMin = q.getMassMin()
            mzMax = q.getMassMax()
            intMin = q.getIntMin()
            intMax = q.getIntMax()
            numVals = q.getNumVals()
            numUsed = q.getNumUsed()

            rt = q.getRetentionTimes()

            if convertToMinutes:
                rt = float(rt) / 60

            mz = resfile.getObservedMass(i)
            charge = q.getCharge()
            ms2mz, ms2int = [], []
            num_peaks = q.getNumberOfPeaks(1)

            for j in range (1, 1+ num_peaks) :
                ms2mz.append(q.getPeakMass(1, j))
                ms2int.append(q.getPeakIntensity(1, j))

            if charge == '2+': # for testing
                query = mgfData(i, title, mzMin, mzMax, intMin, intMax, numVals, numUsed, ms2mz, ms2int, rt, mz, charge)
                query_data.append(query)

    return query_data

def check_residues(peptide, options):
    '''
    Test if > 0 user specified reactive residues is present in a given candidate pep seq
    '''

    if not options.reactiveRes: return False

    reactive_residues = options.reactiveRes
    res_count = 0
    for residue in str(peptide):
        if residue in reactive_residues: res_count += 1
    if res_count > 0:
        return True
    else:
        return False

def filter_peptides(HT_hit, peptides, options):
    '''
    Find mascot peptides that are within defined mz and rt boundaries of a specified HT hit
    '''
    hth_mz = HT_hit.mz
    hth_rt = HT_hit.rt
    charge = HT_hit.charge

    mz_low_limit = float(hth_mz) - options.mzBand[1]/charge #200/2 # TODO - user variables
    mz_high_limit = float(hth_mz) - options.mzBand[0]/charge #70/2 # TODO -divide by 2 for +2 charge state

    filtered_peptides = []

    for p in peptides:
        if p.pep_score > p.identity:
            if float(p.mz) > mz_low_limit and float(p.mz) < mz_high_limit:
                filtered_peptides.append(p)

    return filtered_peptides

def find_mgf_query(HT_hit, mgfdata, HT_ms2_mz_tol, HT_ms2_rt_tol):
    '''
    Find the MGF query associated with the HT hit
    --> this is the MS/MS spectrum of the modified peptide

    This is trivial when HT and DDA data are provided by the same LCMS run
    When separate MS1 and MS2 runs are conducted:
        - Need to correlate HT hits with MS2s of the HT hit accounting for likely RT offsets
            - for hits in different RT ranges

    Also should account for m/z differences induced by MS2 selection of Heavy peptide compared to Light-indexed
    HT hits
    '''

    hth_mz = float(HT_hit.mz)
    hth_rt = float(HT_hit.rt)

    # nested lists of possible queries [mz diff, rt diff, mgf object]
    # queries = []

    # remember hit that minimises mz difference
    min_mz_hit = None
    min_mz_diff = None

    # might make these user variables
    mz_tol = HT_ms2_mz_tol
    rt_tol = HT_ms2_rt_tol

    for mgf in mgfdata:
        mgf_mz = float(mgf.mz)
        mgf_rt = float(mgf.rt)

        mz_diff = math.fabs(hth_mz - mgf_mz)
        rt_diff = math.fabs(hth_rt - mgf_rt)

        if mz_diff < mz_tol and rt_diff < rt_tol:
            #queries.append([mz_diff, rt_diff, mgf])

            if min_mz_hit is not None and mz_diff < min_mz_diff:
                min_mz_hit = mgf
                min_mz_diff = mz_diff
            elif min_mz_hit is None:
                min_mz_hit = mgf
                min_mz_diff = mz_diff

    # find query that minimises mz_diff
    return min_mz_hit

    #return queries # Select first elemnt for testing FIXME

def remove_duplicates(peptides):
    '''
    Remove duplicate entries from Mascot peptide assignment data
        - For each peptide that has been fragmented multiple times, pick entry with highest score
    '''

    rt_tol = 10

    unique_peptide_positions = []

    unique_peptides = []

    for count, i in enumerate(peptides):
        if count % 10 == 0:
            print 'Removing duplicates for peptide %s of %s' %(count, len(peptides))

        imz = i.mz
        irt = float(i.rt)

        pepstats = [i.mz, i.sequence]
        if pepstats not in unique_peptide_positions:
            unique_peptide_positions.append(pepstats)

            duplicates = [i]

            for j in peptides:
                jmz = j.mz
                jrt = float(j.rt)

                rtDelta = math.fabs(irt - jrt)

                if imz == jmz and rtDelta < 10 and irt != jrt:
                    duplicates.append(j)

            # get highest scoring peptide
            duplicates.sort(key = lambda x: float(x.score), reverse = True)
            unique_peptides.append(duplicates[0])

    return unique_peptides

class PeptideFragments(object):
    def __init__(self, residue, residue_index, CRM_mass, a = None, b = None, c = None, x = None, y = None, z = None, correlation_score = None):

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

    def __repr__(self):
        return ( 10 * '%s, ' %(self.residue, self.residue_index, self.CRM_mass, self.correlation_score, self.a, self.b, self.c, self.x, self.y, self.z))

def get_fragment_mod_masses(site, frag_mod_str, var_mods, ht_mz, pep_mz):
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
    pep_mz += assigned_mod_mass / 2

    # CRM mass
    mod_mass = (float(ht_mz) - float(pep_mz)) *2

    # get masses of mods referenced in mod_string
    frag_mod_mass = []
    for i, resMod in enumerate(frag_mod_str):

        # if residue is not the candidate CRM site - check for native mods
        if i != site:
            # if value in frag_mod_str == 0 ---> No modofication present at this site in Mascot peptide
            if resMod == '0': frag_mod_mass.append(float(0))

            # if valeu != 0 ---> Some Mascot-assigned modification is present - find mass of this mod in var_mods list
            else:
                for mod in var_mods:
                    if str(mod.modIndex) == resMod:
                        frag_mod_mass.append(mod.modDelta)

        # if residue is candidate CRM site - check for native mods - if present, adjust apparrent CRM mass accordingly
        elif i == site:

            if resMod == '0':
                # no native modification assigned - append full mod mass
                frag_mod_mass.append(mod_mass)

            else:
                # native mod present --> (ht_mz - pep_mz) * charge represents the NET CHANGE in mass
                # add mass of native mod to observed ht-pep mass difference to get full CRM candidate mass

                # observed difference between HT and Mascot peptides
                ht_pep_diff = (float(ht_mz) - float(pep_mz)) *2

                # mass of native mod in mascot peptide
                native_mod_mass = None
                for mod in var_mods:
                    if str(mod.modIndex) == resMod:
                        native_mod_mass = mod.modDelta

                mod_mass = ht_pep_diff + native_mod_mass
                frag_mod_mass.append(mod_mass)

    return frag_mod_mass, mod_mass

def rolling_modification_fragments (peptide, mod_string, var_mods, ht_mz, pep_mz, reactiveRes, types = ('y'), maxcharge = 2):
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

    frag_objects = []

    # move CRM modification site along peptide
    for site in range(len(peptide)):

        if peptide[site] not in reactiveRes: continue

        # get the mod string for residues in this fragment
        frag_mod_str = list(mod_string)

        # calculate mass of unmodified peptide
        calc_pep_mz = mass.fast_mass(peptide, charge = 2)

        # create a list of masses to add/subtract from each residue
        frag_mod_mass, CRM_mass = get_fragment_mod_masses(site, frag_mod_str, var_mods, ht_mz, calc_pep_mz)

        pepFrags = PeptideFragments(peptide[site], site, CRM_mass)

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

                        b.append([peptide[:i], float(mz), ion_type])

                    if ion_type in 'xyz':

                        # generate pure sequence of fragment
                        frag = peptide[i:]

                        # get mass of relevant mods
                        mods = frag_mod_mass[i:]

                        # calculate mass of base fragment
                        mz = mass.fast_mass(peptide[i:], ion_type = ion_type, charge = charge)

                        # add total mass of modifications - including CRM
                        mz = mz + sum(mods)

                        y.append([peptide[i:], float(mz), ion_type])


        if len(a) != 0: pepFrags.a = a
        if len(b) != 0: pepFrags.b = b
        if len(c) != 0: pepFrags.c = c
        if len(x) != 0: pepFrags.x = x
        if len(y) != 0: pepFrags.y = y
        if len(z) != 0: pepFrags.z = z
        frag_objects.append(pepFrags)

    if len(frag_objects) == 0:
        return None
    else:
        return frag_objects

def similarity_score(expt_mz, expt_int, theoretical_peptide_set, options):
    '''
    Determine the correlation between the theoretical ion series
    (derived form mascot-assigned sequence) and the MS/MS spectrum of a HT hit
    '''

    maxScore = 0
    maxScorePeptide = None

    # normalise experimental intnsities
    maxInt = float(max(expt_int))
    expt_int = [x/maxInt * 100 for x in expt_int]

    results = []

    if options.matchTol:
        tol = options.matchTol
    else:
        tol = 0.2

    #     #ofx = open('fit_data.dat','a')
    for m, theorPeptide in enumerate(theoretical_peptide_set):

        frags = [theorPeptide.a, theorPeptide.b, theorPeptide.c, theorPeptide.x, theorPeptide.y, theorPeptide.z]
        frags = [x for x in frags if x is not None]

        # flatten mz array
        #theor_mz = [i[1] for l in frags for i in l]
        theor_mz = [i for l in frags for i in l]
        explained_signal = 0
        num_matched = 0

        ionTracker = {
            'a': 0,
            'b': 0,
            'c': 0,
            'x': 0,
            'y': 0,
            'z': 0
            }

        for e_index, e_mz in enumerate(expt_mz):
            for t_index, t_mz in enumerate(theor_mz):
                if float(e_mz) < float(t_mz[1]) + tol and float(e_mz) > float(t_mz[1]) - tol:
                    explained_signal += expt_int[e_index]
                    num_matched += 1

                #    if ionTracker[t_mz[2]] == t_index - 1:
                #        ionTracker[t_mz[2]] += 1

        #multiplier = sum(ionTracker.values())
        score = explained_signal * num_matched*2#multiplier

        # add correlation_score attribute

        theorPeptide.correlation_score = score
        results.append(theorPeptide)

    return results

class formulas(object):
    def __init__(self, formula, mass, ppm, formula_dict = None):
        self.formula = formula
        self.mass = mass
        self.ppm = ppm

        if formula_dict:
            self.formula_dict = formula_dict
        else:
            self.formula_dict = self.get_formula(formula)

    def get_formula(self, f_string):
        '''
        Create a dictionary holding the molecular formula
            - keys = elements
            - values = stoichiometry
        '''

        f_string = f_string.replace('*','')
        i = re.findall(r'([A-Z][a-z]*)(\d*)', f_string)

        formula = {}
        for x in i:
            if x[1] != '' and x[1] != '0':
                formula[x[0]] = x[1]

        formula = formula
        return formula

    def __repr__(self):
        return ( 3*'%s ' %(self.formula, self.mass, self.ppm))

def make_formula_str_from_dict(fdict):

    mf_string = ''
    for k,v in fdict.iteritems():
        mf_string += '%s%s ' %(k,v)

    return mf_string.strip() # remove trailing white space

def get_possible_formula_masses(options):
    '''
    Calculate all possible molecular weights using allowed atom ranges and store in numpy array

    Return Value:
        1) NP array of formula mass values
        2) list of dicts corresponding to NP array indices
    '''

    # data structure: values stored in numpy array, corresponding formulae stored in list
    # --- NB: np arrays with mixed datatypes are possible - this might be faster - TODO
    data = []
    data_list = []

    # append function
    add = data_list.append

    atoms, ranges = [],[]
    for k,v in options.atomDict.iteritems():
        atoms.append(k)
        ranges.append(v)

    combinations = itertools.product(*ranges)

    # define datatype to be used in numpy arrays
    datatype = 'float, str'

    # calculate mass of each combination
    for x in combinations:
        mf = dict()

        for i, atom in enumerate(atoms):
            mf[atom] = x[i]

        mz = mass.calculate_mass(composition = mf)

        # add value to array and formula to list
        data.append(float(mz))
        add(mf)

    # convert data list to numpy array
    data = np.asarray(data)

    return data, data_list

def find_formula_in_formula_list(formula_array, formula_list, target, options):
    '''
    Find formulae that are within the tolerances of a given target mass in the formula array

    Return value: list of formula objects:
        - Attributes:
            MF = molecular formula string
            formula_dict = MF in dictionary form: keys = elements, values = stoichiometry
            mz = mass
            ppm = ppm error between formula and target
    '''

    # calculate mass tolerance
    #tol = 0.5    # TODO user variable

    tol = float(options.ppmTol) * target / 1000000

    # calculate boundaries used to search formula mass array
    lower_boundary = target - tol
    upper_boundary = target + tol

    # list of 'formula' objects
    matches = []

    # find array elements that fit within boundaries - returns a touple of arrays --> select first element
    match_indices = np.where((formula_array > lower_boundary) & (formula_array < upper_boundary))[0]

    # instantiate formula objects and add to matches list
    for match in match_indices:
        MF_dict = formula_list[match]
        mz = formula_array[match]
        ppm = (float(target) - float(mz))/float(mz) * math.pow(10, 6)
        if math.fabs(ppm) < 100:
            matches.append(formulas( make_formula_str_from_dict(MF_dict), mz, ppm, formula_dict = MF_dict))

    # sort matches list by ppm error
    matches.sort(key = lambda x: math.fabs(float(x.ppm)))
    return matches

def get_residual_mass_error(fragment_formulae, CRM_formulae, maxRME):

    '''
    Calculate RME for hits of a given formula for a given HT/mascot peptide pair

    fragment_formulae > list of cStructure elements
                      > attributes are: self.structure_type, self.smile, self.MW, self.formula_plain, self.formula_dict

    CRM_formulae > list of formulas elements (defined in NTPS)
                 > attributes are: self.formula, self.mass, self.ppm

    Two attributes are added to CRM_formula object:
        - RME - value of the minimal RME for the given CRM
        - RME_frag - cStructure object defining the fragment with the minimum RME

    '''
    returnFormulae =  []

    # iterate through possible CRM formulae
    for CRM in CRM_formulae:

        RMS = None
        RMS_frag = None

        # Calculate RME between the CRM and each possible fragment
        for frag in fragment_formulae:

            residual_mass_error = 0

            # get formula dictionaries - Keys = atoms, Values = Stoichiometry
            CRM_dict = CRM.formula_dict
            frag_dict = frag.formula_dict

            # get list of atoms in frag and CRM
            CRM_keys = set(CRM_dict.keys())
            frag_keys = set(frag_dict.keys())

            # elements that occur at least once in either CRM_keys or frag_keys sets
            all_elements = CRM_keys.union(frag_keys)

            # calculate differences between atoms that are common to both CRM and frag
            # TODO --> should non-common elements be included as well?
            for element in all_elements:

                # if element in both lists - calc difference and add mass to RME
                if element in CRM_keys and element in frag_keys:
                    element_error = int(math.fabs(int(CRM_dict[element]) - int(frag_dict[element])))
                    if element_error > 0:
                            # calculate mass and add to RME
                            mz = mass.calculate_mass(formula = '%s%s'%(element, element_error))
                            residual_mass_error += mz

                # else add mass to RME
                elif element in frag_keys:
                    element_error = int(frag_dict[element])
                    mz = mass.calculate_mass(formula = '%s%s'%(element, element_error))
                    residual_mass_error += mz

                elif element in CRM_keys:
                    element_error = int(CRM_dict[element])
                    mz = mass.calculate_mass(formula = '%s%s'%(element, element_error))
                    residual_mass_error += mz

            if RMS is not None:
                if residual_mass_error < RMS:
                    RMS = residual_mass_error
                    RMS_frag = frag
            else:
                RMS = residual_mass_error
                RMS_frag = frag

        if RMS < maxRME:
            CRM.RME = RMS
            CRM.RME_frag = RMS_frag.img
            returnFormulae.append(CRM)

    return returnFormulae

def generate_fragemnts(drugSMILES):
    '''
    Disconnect rotatable bond of SMILES specified molecule
    Save fragments calculate formulae
    Load fragments into structure viewer widget
    '''

    # Clear class level list
    molecule = []

    # Append cStructure mol to class level list
    molecule.append(cStructure('precursor', drugSMILES, None))

    # SMARTS string for bond disconnection
    patt = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]')

    # Init parent mol object
    mol = Chem.MolFromSmiles(drugSMILES)

    # find the rotatable bonds
    bonds = mol.GetSubstructMatches(patt)

    # create an editable molecule, break the bonds, and add dummies:
    #all_smis = [drugSMILES]

    # disconnect rotatable bonds
    for a,b in bonds:
        em = Chem.EditableMol(mol)
        nAts = mol.GetNumAtoms()
        em.RemoveBond(a,b)
        em.AddAtom(Chem.Atom(0))
        em.AddBond(a,nAts,Chem.BondType.SINGLE)
        em.AddAtom(Chem.Atom(0))
        em.AddBond(b,nAts+1,Chem.BondType.SINGLE)
        nAts+=2
        p = em.GetMol()
        Chem.SanitizeMol(p)
        smis = [Chem.MolToSmiles(x,True) for x in Chem.GetMolFrags(p,asMols=True)]
        for smi in smis:
            #all_smis.append(smi)
            molecule.append(cStructure(None, smi, None))

    return molecule

def write_output(results, options):

    of1 = open(options.outFile, 'wt')
    of1.write('# %s: %s\n' %('htInput', options.mascotInput))
    of1.write('# %s: %s\n\n' %('mascotInput', options.mascotInput))
    of1.write('# index, pep_sequence, mod_mass, pep_mz, pep_rt, ht_mz, ht_rt, score, modified_residue, modified_residue_index\n')
    for i, r in enumerate(results):
        of1.write('<BEGIN>\n')
        of1.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n' %(i, r.pep_sequence, r.mod_mass, r.pep_mz, r.pep_rt, r.ht_hit_mz, r.ht_hit_rt, r.score, r.modified_residue, r.modified_residue_index))
        of1.write('<FORMULAE>\n',)
        for f in r.formulae:
            of1.write('%s, %s, %s\n'%(f.formula, f.ppm, f.RME))
        of1.write('<END>\n\n')

    of1.close()
    return

def main(options, q = None):
    '''
    TODO - need to include a print function for when NTPS is not run through GUI
    '''
    startTime = datetime.now()

    HT_data = []

    if options.drugSMILES:
        molecule = generate_fragemnts(options.drugSMILES)
    else:
        molecule = options.molFrags

    # iteratively parse HT files
    for htInput in options.htInput:

        htInputFile = htInput[0] # target file name
        htInputCharge = htInput[1] # target file charge

        additional_HT_attributes = dict()
        additional_HT_attributes['charge'] = int(htInputCharge)

        headers, this_HT_data = read_hitime_files(
                                            htInputFile,
                                            addHtAttributes = additional_HT_attributes
                                            )

        # check if this ht data is a verified peak list - default to false if not
        peakList = headers.get('validated', False)

        if not peakList:
            # run automated peak picking if not a peak list
            this_HT_data = get_HT_regions(this_HT_data,
                                    0.3, # Drt
                                    0.1, # Dmz
                                    minScore = options.minHTScore
                                    )
        HT_data += this_HT_data



    print('%s HiTIME hits found' %len(HT_data))

    ''' Load MS/MS and Mascot data '''
    # Load MS/MS query data
    mgfdata = get_mgf_data(options.mascotInput, convertToMinutes = options.toMins)
    print('%s MGF queries found' %len(mgfdata))

    # get mascot-assigned peptides and variable mods used in fitting
    # 'peptides' is a lsit of mascotHit objects
    # 'var_mods' is a list of varMods objects

    print 'Reading Mascot input file'
    resfile = load_mascot_results(options.mascotInput[0])

    print 'Reading peptide results'
    peptides, var_mods = get_peptide_results(resfile, convertToMinutes = options.toMins)

    # remove duplicates in peptide assignment data
    peptides = remove_duplicates(peptides)

    print('%s peptides found' %len(peptides))
    print('%s var_mods found' %len(var_mods))

    '''
    peptide objects have an attribute 'mod_string'
        var_mod objects have the attributes: modIndex, modName, modDelta, modNeutralLoss
        in mod string - 0 indicated no modification. Otherwise the modification is specified by an integer corresponding to the modIndex
        typically, modNeutralLoss == 0 as no atoms are lost upon modifications
    '''

    # caluclate possible CRM masses from allowed formulae
    '''
    TODO:
        need to pass atom range dictionary to this function
    '''
    formula_array, formula_list = get_possible_formula_masses(options)

    count = 1
    results = []
    found_sequences = 0

    # iterate through HT hits
    for index, ht_hit in enumerate(HT_data):

        if index % 2 == 0:
            print 'Processing hit %s of %s' %(index, len(HT_data))

        # create results list for each HT hit -- append all results then return the highest scoring hit
        HT_hit_results = []

        # filter assigned peptides based on mz and rt tolerances
        peptide_subset = filter_peptides(ht_hit, peptides, options)
        sequences_done = []

        # check if any assigned peptides are found
        if len(peptide_subset) == 0 : continue

        # find ms2 queries near HT hits

        nearest_mgf_query = find_mgf_query(
                                            ht_hit,
                                            mgfdata,
                                            options.HT_ms2_mz_tol,
                                            options.HT_ms2_rt_tol
                                            )

        if nearest_mgf_query is None:
            continue # no MS2 spectrum correlating to the HT hit

        ht_hit.mz = nearest_mgf_query.mz

        max_score = 0 #

        #of1 = open('sequences.dat','w')
        for counter, peptide in enumerate(peptide_subset):

            # test that user defined reactive residues are in peptide
            # if check_residues(peptide.sequence) is True:

            if not check_residues(peptide.sequence, options): continue

            # get sequence and determine theoretical fragments
            theoretical_peptide_set = rolling_modification_fragments(peptide.sequence, peptide.mod_string, var_mods, ht_hit.mz, peptide.mz, options.reactiveRes, types = options.ionTypes)

            expt_mz = nearest_mgf_query.ms2mz
            expt_int = nearest_mgf_query.ms2int
            expt_int = [ float(x) for x in expt_int]

            # calculate correlation score between spectra
            # find theoretical peptide with highest score
            result = similarity_score(expt_mz, expt_int, theoretical_peptide_set, options)

            if result == None: continue
            else:
                count += 1
                result.sort(key = lambda x: int(x.correlation_score), reverse = True)

                result = result[:1] #result[:1]

                for theor_peptide in result:

                    score = theor_peptide.correlation_score
                    CRMMass = theor_peptide.CRM_mass

                    # get list all fragment ions produced by rolling mods function
                    theor_mz = []
                    for i in [theor_peptide.a, theor_peptide.b, theor_peptide.c, theor_peptide.x, theor_peptide.y, theor_peptide.z]:
                        if i is not None: theor_mz = theor_mz + i

                    theor_mz = [x[1] for x in theor_mz]
                    theor_int = [nearest_mgf_query.intMax*1.2 for x in theor_mz]

                    # get formulae that match the mz difference in the peptide pair
                    formulas = find_formula_in_formula_list(formula_array, formula_list, CRMMass, options)

                    # calculate RME and apply RME cutoff
                    formulas = get_residual_mass_error(molecule, formulas, options.maxRME)

                    if len(formulas) != 0 and int(score) > 0:

                        delta_rt = float(ht_hit.rt) - float(peptide.rt)
                        result = Match(peptide.sequence, peptide.mz, peptide.mod_string, ht_hit.mz, CRMMass, score, formulas, count, np.asarray(theor_mz), np.asarray(theor_int), np.asarray(expt_mz), np.asarray(expt_int), peptide.identity, peptide.homology, peptide.pep_score, ht_hit.score, theor_peptide.residue, theor_peptide.residue_index, delta_rt, ht_hit.rt, peptide.rt)

                        # append result HT_hit_results
                        HT_hit_results.append(result)

                        #results.append(result)
                        found_sequences += 1


        # sort hits by score and return to guiprocs
        if len(HT_hit_results) != 0:
            #done = []
            HT_hit_results.sort(key = lambda x: float(x.score), reverse = True)
            scores = [x.score for x in HT_hit_results]

            this_result = HT_hit_results[0]
            this_result.corrs = scores
            try:
                this_result.second = HT_hit_results[1]
            except IndexError:
                pass

            try:
                this_result.third = HT_hit_results[2]
            except IndexError:
                pass

            try:
                this_result.third = HT_hit_results[4]
            except IndexError:
                pass

            results.append(this_result)

            # for counter, res in enumerate(HT_hit_results):
            #     if res.mod_mass not in done:
            #         if counter > 0: continue
            #         done.append(res.mod_mass)
            #         res.corrs = scores
            #         results.append(res)

    print 'Execution time is: %s' %(datetime.now() - startTime)

    print 'Number of correlations is: ', len(results)

    results.sort(key = lambda x: float(x.score), reverse = True)

    if options.pickle: import pickle

    if options.outFile:
        write_output(results, options)

        # for testing
        if options.pickle:
            with open(r'%s' %(options.outFile + '.pickle'), 'wb') as f:
                pickle.dump(['done', results, HT_data], f)

    returnResults = []
    for i, r in enumerate(results):
        r.count = i
        returnResults.append(r)

    if q is not None:
        q.put(['done', results, HT_data])

    else:
        return returnResults

class Match(object):
    def __init__(self, pep_sequence, pep_mz, pep_mod_str, ht_hit_mz, mod_mass, score, formulae, index, theor_mz, theor_int, expt_mz, expt_int, identity, homology, mascot_score, HT_score, modified_residue, modified_residue_index, delta_rt, ht_hit_rt, pep_rt):
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
        self.identity = identity
        self.homology = homology
        self.mascot_score = mascot_score
        self.HT_score = HT_score
        self.modified_residue = modified_residue
        self.modified_residue_index = modified_residue_index
        self.delta_rt = delta_rt
        self.ht_hit_rt = ht_hit_rt
        self.pep_rt = pep_rt

    def __repr__(self):
        return ( 14 * '%s,' %(self.pep_sequence, self.pep_mz, self.pep_mod_str, self.ht_hit_mz, self.mod_mass, self.score, self.formulae, self.index, self.theor_mz, self.theor_int, self.expt_mz, self.expt_int, self.identity, self.homology, self.mascot_score, self.HT_score, self.modified_residue, self.modified_residue_index, self.delta_rt, self.ht_hit_rt, self.pep_rt))

class GUIOptions(object):

    def __init__(self, args):
        for k, v in args.iteritems():
            setattr(self, k, v)

def apply_defaults(args):
    return

def params(q, args): #args):
    '''
    Init function called when script is executed from within the GUI
    '''

    options = GUIOptions(args)
    results = main(options, q)
    return results

if __name__ == '__main__':
    options = parser.parse_args()

    options.htInput = [options.htInput]
    options.mascotInput = [options.mascotInput]

    guiMode = False
    sys.exit(main(options))
