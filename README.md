# Introduction
This software provides various data analysis tools that can be used to interpret combined twin-ion / shotgun proteomics experiments. The original aim of these experiments has been to identify the protein targets of reactive drug metabolites formed during liver microsome incubations.

The software utilised the High-resolution twin-ion metabolite extraction (HiTIME) algorithms developed by Dr. Bernard Pope and Dr. Andrew Isaac at the University of Melbourne. A publication describing implementation of HiTIME can be found <a href = 'http://pubs.acs.org/doi/abs/10.1021/ac504767d'> here </a>.  

# Installation
Xenophile was developed and deployed on Ubuntu linux 14.04. Installation instructions a fresh install of Ubuntu are provided below:

```bash

# install foundation stuff
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install python-pip
sudo pip install --upgrade pip
sudo apt-get install python-setuptools python-dev build-essential automake autoconf libtool git

# install scipy stack
sudo apt-get install python-tk
sudo pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose

# install pyteomics and pymzml
sudo pip install lxml image pyteomics 

# install RDkit
sudo apt-get install python-rdkit librdkit1 rdkit-data

# install libspatialindex rtree
git clone https://github.com/libspatialindex/libspatialindex.git
cd libspatialindex/
./autogen.sh
./configure; make; make install
sudo apt-get install libspatialindex-dev
ldconfig
sudo pip install rtree
cd

# install qt4
sudo apt-get install  libglew-dev libcheese7 libcheese-gtk23 libclutter-gst-2.0-0 libcogl15 libclutter-gtk-1.0-0 libclutter-1.0-0

# install pyqt4
sudo apt-get install python-qt4 qt4-dev-tools python-qt4-dev pyqt4-dev-tools

# install xenophile
git clone  https://github.com/mgleeming/Xenophile.git
cd Xenophile
sudo python setup.py install
cd

```

**HiTIME**  
   The HiTIME algorithm performs the initial LC/MS data weighting by correlation with the twin-ion signature. Two iterations of this software are available.

- **Initial python 2.7 release** - This is distributed with the Protein-CRM identification code. Provides basic HiTIME functionality but is relatively slow. Documentation is available <a href = 'https://github.com/bjpop/HiTIME'> here </a>.
- **Production C++ version** </a> - Significantly faster with enhanced functionality and OpenMS integration. Requires OpenMS which must be built from source. HiTIME-CPP documentation and OpenMS build instructions for Mac OSX and Linux are available <a href = 'https://github.com/bjpop/HiTIME-CPP/tree/local_maxima'>here </a>. Xenophile will use HiTIME-CPP by default if available.  

**Python packages**
- **pymzml** - Reading mzML mass spectrometry data files
- **rtree** - Spatial indexing and peak-picking
- **numpy** - Array math
- **rdkit** - Manipulation of chemical structures
- **pyteomics** - Calculation of molecular masses of peptides and chemical formulae
- **PyQt4** - GUI presentation and functions
- **pyqtgraph** - Plotting library optimised for GUI applications  


**Other**
- **mascotparser** - Parsing of raw Mascot '.dat' output files. Note that a version of mascotparser compatible with linux operating systems is distributed with the xenophile source. For use on non-linux systems, the appropriate mascotparser distribution can be obtained free-of-charge<a href = 'http://www.matrixscience.com/msparser.html'> here</a>. Ensure that the python directory within the mascotparser source is added to your system PATH prior to use.


# HiTIME search
Performs initial LC-MS data weighting by correlation with the twin-ion signature. This is the starting point for protein-CRM adduct detection. For best results, it is recommended that HiTIME searching be performed on mzML files that contain regularly spaced MS1 spectra - i.e. mzML files for data-dependent acquisition LC-MS/MS runs result in significantly lower scores for true twin-ion signals due to irregular timing between MS1 acquisitions.

Modules included in this section are:
**HiTIME search** - Define input/output files, parameters and execute search.

**Background subtraction** (optional) - point-by-point subtraction of a background or contol data file from a treatment data file.

  * **__m/z__ Tolerance** - m/z window used to align HiTIME peaks for subtraction.  
  * **RT Tolerance** - Retention time window used to align HiTIME peaks for subtraction. Time units are those present in the HiTIME data file.
  * **Score Cutoff** - Minimum HiTIME score required for data points to be considered for subtraction.  

**Results post-processing** - Facilitates analysis of HiTIME results by retrieving Extracted Ion Chromatograms (EICs) and mass spectra for each candidate peak.
  * **Minimum HT Score** - Minimum HiTIME score required for data points to be considered  
  * **m/z Width** - m/z window used to align HiTIME peaks for subtraction.  
  * **Retention Time Width** - Retention time window used to align HiTIME peaks for subtraction.  
  * **RT exclusion** - Optional. Exclude duplicate hits at longer retention times. Prevents multiple extractions of the same compound due to sample bleed down the chromatographic column.  
  * **mzDelta** - m/z difference between heavy and light compounds.  
  * **EIC width** - m/z range used to produce extracted ion chromatogram for identified HiTIME hits (+/- m/z)  

**Results Visualisation** - Framework to facilitate visual inspection of the EIC traces and mass spectra obtained from data postprocessing.

# Non-targeted CRM identification

The non-targeted peptide identification algorithm is designed to elucidate the formula of a CRM *de novo*.

#### Input structure specification
**SMILES String**: SMILES string defining the structure of the administered drug molecule. The string can be pasted into the input field or, alternatively, the molecule can be entered using a graphical structure editor that can be accessed by pressing the 'draw structure' button.

  After defining the input molecule, pressing the 'Generate Fragments from Structure' button will recursively cleave the molecule at each rotatable bond. The starting structure, and the fragments generated by this process, can be reviewed in the adjacent 'Structure Viewer' panel.

**Element Ranges**: Specify the allowed element ranges used for determining CRM stoichiometries. These determine the maximum and minimum number of each element that is allowed to be used in composing candidate CRM formulae. A combination of discrete value(s) and ranges can be specified as follows:

| Element | Range specification  | Result     |
| ------- | ------------- | -------------     |
| C       | 6-10          | 6, 7, 8, 9, 10    |
| C       | 3, 6-10       | 3, 6, 7, 8, 9, 10 |
| C       | 3, 5, 7       | 3, 5, 7           |
| C[13]   | 0, 6          | 0, 6              |
| O[18]   | 0, 2          | 0, 2              |

If not specified otherwise, the most abundant naturally occurring isotope is used in mass calculations. To include other isotopes, enclose the isotope in square brackets following the atomic symbol.

#### Input data specification

**Mascot Threshold** - Minimum peptide confidence interval to be used in HiTIME correlation. (identity > homology)

**Charge State Range** - Range of assigned peptide charge states to be considered for HiTIME correlation. For tryptic digests, +2 to +4 is recommended.

**HT-MS2 correlation**:
HiTIME hits identified from MS1 data need to be linked to associated MS2 spectra. The MS2 spectra are extracted from peak lists retained in the Mascot .dat file. This is trivial when the MS1 and MS2 data in use comes from the same LC/MS run as these data points will overlap. However when MS1 data (used for HT scoring) and MS2 data are collected in different LC/MS runs, retention time offsets between HT peaks and associated MS2 spectra need to be taken into account.

To account for these potential offsets, an extra parameter set is included.

- **HT MS2 m/z offset** - mz tolerance between HiTIME peaks and associated MS2 spectrum
- **HT MS2 RT offset** - rt tolerance between HiTIME peaks and associated MS2 spectra

    NOTE: the retention time offset is the most important here. If MS1/MS2 are in different runs, make this number large enough to capture the HT hit and the associated MS2. The software will select mascot peptides that minimise the RT difference between HT and MS2 peaks.  


#### MS/MS correlation parameters

**Match Ion Types**: Peptide sequence ion types that should be considered when correlating MS2 spectra of HiTIME hit peptides and native peptides. In general, 'b' and 'y' type ions should be selected for thermal fragmentation methods and 'c' and 'z' should be used for fragmentation via ETD or ECD.

**Match Tolerance**: m/z error (Da) allowed when correlation MS2 spectra of HiTIME and native peptides

**Reactive Residues**: Peptide residues that are to be considered 'reactive' toward modification with potential CRMs. These selections can be used to restrict the search space in MS2 correlation. If none are selected, all residues will be considered reactive.

#### CRM Parameters

**_m/z_ band**: _m/z_ range relative to a HiTIME hit that is searched for potential native counterpart peptides. Initial values for this parameter are derived from the elemental compositions specified above.

  For example, _m/z_ band specification of "75-200" will result in the following; if a HiTIME hit is detected at _m/z_ = 500, the _m/z_ range searched for matching native peptides will be _m/z_ 300 to _m/z_ 425 which is equivalent to 500 - 300 to 500 - 75.

**RT band**: optional parameter. If specified, only falling within this range will be considered as potential native counterparts.

**ppm tolerance**: Return all possible formulae that match the candidate CRM mass to within the specified ppm tolerance.

# Targeted protein-CRM adduct search

Correlates MS/MS data to find peptides that are observed as twin-ions and are also assigned as CRM adducts by mascot searching.

#### Input data specifications

**Mascot**: Mascot search output file.  
**HiTIME**: HiTIME search output file. These can be either raw output files or validated HiTIME hit lists. Multiple HiTIME data files can be added to allow for explicit specification of target lists for multiple charge states. The charge state corresponding to each provided HiTIME file must be specified in the 'charge' column of the HiTIME data file table.

#### CRM specification

All fixed and variable modifications found in the mascot output file provided in the 'input data specification' field are presented in the 'All modifications' table. Select the modifications that correspond to the CRMs used in the mascot search and press 'add' to move these to the 'CRM modifications' list.

The correlation algorithm will search for peptides that have these modifications and are also observed as twin-ion hits.
