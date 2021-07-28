## PIPPy
## A Python and Fortran code for automatically fitting permutationally invariant polynomials

## Authors
Daniel R. Moberg and Ahren W. Jasper\
Argonne National Laboratory, 2021\
Web: github.com/Auto-Mech/PIPPy\
Email: ajasper@anl.gov

References:\
(1) D. R. Moberg and A. W. Jasper, PIPPy, Argonne National Laboratory, 2021.\
(2) D. R. Moberg and A. W. Jasper, Permutationally invariant polynomial expansions with unrestricted complexity, J. Chem. Theory Comput., submitted (2021)

---------------------------------------------------------------------------------------
## I. DISTRIBUTION AND INSTALLATION

PIPPy is distributed with four subdirectories:\
- doc/: Contains this manual in txt format\
- python/: Contains the Python wrapper\
- runs/: Contains several example input and output files\
- src/: Contains the Fortran code and a compilation script named compile.sh

To install PIPPy, first compile the Fortran executable by running ./compile.sh in the 
src/ directory. The compiled executable is named pip.x and appears in the src/ directory.

The Fortran code can be run independently of the Python wrapper with a properly prepared
input file. An example of running a parallel executable is:\
   `mpiexec -n 36 src/pip.x < input > output`

Alternatively, the Python wrapper can be used to set up and run PIPPy. The Python script
requires the Mako Templates module (see [https://www.makotemplates.org/]). The user then
sets up a key/pair list using the keywords defined for the standard input in II.A and
runs (ideally with Python 3.0 or later):\
   `python pip.py input.py`

The Python script will generate the required Fortran input file (named input), run the
code, and collect the standard output in the file named "output".

---------------------------------------------------------------------------------------
## II. INPUT FILES

PIPPy requires a standard input file and at least one data file. All input files should appear in the runtime directory. 

### II.A. STANDARD INPUT

The standard input file contains a series of records, each a single line of input parameters. Below is a list of each record and associated parameters, with data type in parentheses.\
Key: dp = double precision, int = integer, char * _n_ = character * _n_, log = logical


- **Record 1**: DataTrain  DataTest [Example: train.dat test.dat]
    - DataTrain (char * 10): Filename of the in sample (training) data file
    - DataTest (char * 10): Filename of the out of sample (test) data file

    **Note:** The required format for the data files is given in Section II.B.\
           If DataTest is set to "none", no test set evaluation is made.

- **Record 2**: NumWrite  Units [Example: 2 10 12]
    - NumWrite (int): Number of extra output files to write\
        = -1, no extra output files are written.\
	= 0, all output files are written.\
	&gt; 0, read a list of NumWrite extra output files.
    - Units[NumWrite] (int): List of additional output files to write.\
	 1 : Write additional debug information to the standard output\
	10 : (basisinfo.dat) Write detailed basis set information\
	11 : (vtrain.dat) Compare ab intio training energies and fitted energies\
	12 : (vtest.dat) Compare ab intio test energies and fitted energies\
	20 : (xmat.dat) Write the X matrix containing the symmetrized basis functions (columns) evaluated for the training data (rows)

- **Record 3**: RangeParameter  RefEnergy [Example: 100. 15.]
    - RangeParameter (dp): Range parameter used in the training and test set weight function
    - RefEnergy (dp): Reference energy used in the training and test set weight function

    **Note:** Reasonable choices for RefEnergy include the binding energy for van der Waals systems and the saddle point energy for reactive systems, while RangeParameter describes the energy range required for the intended application.

- **Record 4**: NumRanges  EnergyRanges [Example: 3 300. 150. 50.]
    - NumRange (int): Number of energy limits to consider when reporting averaged errors.
    - EnergyRange[NumRange] (dp): Values of energy range limits.

- **Record 5**: NumAtoms [Example: 6]
    - NumAtoms (int): Number of atoms

- **Record 6**: Symbols [Example: C H H H H  H]
    - Symbols[NumAtoms] (char*2): Atomic symbol for each atom

- **Record 7**: AtomGroups [Example: 1 2 2 2 2  2]
    - AtomGroups[NumAtoms] (int): Atomic permutation group for each atom

- **Record 8**: ReadBasis  FactorOrder  TotalOrder  IMode [Example: F 6 6 1]
    - ReadBasis (logical): Set to TRUE to read the basis from basis.dat. Otherwise the basis is generated and written to basis.dat.
    - FactorOrder (int): Maximum allowed order for a single factor
    - TotalOrder (int): Maximum total order allowed for a term
    - IMode (int): Controls various choices for generating the PIP expansion\
        -1 : Use intermolecular distances only when generating terms\
        0 : Use all PIP terms\
        1 : Remove unconnected terms\
        2 : Remove unconnected and intramolecular-only terms\
        3 : Remove intramolecular-only terms

        **Note:** When IMode != 0, one or more fragment groups must be assigned in Records 9 and 10.

- **Record 9**: NumChannels [Example: 1]
    - NumChannels (int): Number of fragment channels

- **Record 10**: FragmentGroups [Example: 1 1 1 1 2  2]
    - FragmentGroups[NumAtoms] (int): Fragment group assignments

An example of the standard input for CH4 + H &lt;=&gt; CH3 + H2:
```
 train.dat  test.dat           ! Training (in-sample) and test (out-of-sample) data sets
 2 10 12                       ! Write two extra output files: basisinfo.dat and vtest.dat
 120. 15.                      ! Range parameter & reference energy used for weigting the data sets
 3 300. 150. 50.               ! Report weighted RMS errors for the energy ranges 300-150, 150-50, and > 50 kcal/mol
 6                             ! Number of atoms
 C H H H H  H                  ! Atom labels
 1 2 2 2 2  2                  ! Exchange all H atoms
 F 6 6 1                       ! Generate a PIP66 basis and exclude unconnected terms
 1                             ! Define one molecular fragment channel
 1 1 1 1 2  2                  ! Remove unconnected terms for CH3 + H2
```

An example of the standard input for parametrizing the intermolecular energy of HO2 + N2:
```
 ai.all ai.test                ! Training (in-sample) and test (out-of-sample) data sets
 -1                            ! Print only standard output files
 650. -650.                    ! Range parameter & reference energy used for weigting the data sets
 4 12000. 4000. 0. -100.       ! Energy ranges for error analysis (cm-1, should match on training and test data)
 5                             ! Number of atoms
 O O H  N N                    ! Atom labels
 1 1 2  3 3                    ! Atomic permutation group assignments
 F 3 3 1                       ! Read basis flag, xi(factor), xi(total), IMode
 1                             ! Number of fragment channels
 1 1 1  2 2                    ! Fragment channel assignment
```

### II.B. TRAINING AND TEST SETS

These data files contain the training and test sets. PIPPy does not currently include functionality for creating these files. Both files follow the same format:

Record 1:                   nconfig ! Number of data\
Record 2:                   natom   ! Number of atoms\
Record 3:                   iconfig   dummy   vai  ! Index, dummy(not used), ab initio energy\
Record 4 to 4+natom-1:      symb   x   y   z ! Atomic symbol, Cartesian coordinates\
Repeat records 2 to 4+natom-1 nconfigs times.

**Note:** Units are never converted such that the generated expansions should provide energies and forces in the same units as the training data. One caveat: The Morse range parameter is hard coded as unity in whatever units are used for the geometry in the training data. This choice is more appropriate when using Angstroms than when using other units for distance.

An example of a training data file for CH4 + H:
```
 252889
           6
           1   10.560402393341064        53.686826233325540
 C    1.5173673060615555E-002   4.4702465630787759E-002  -1.4358965674743889E-002
 H   0.91590037301702620       0.90465104728081480      -0.20951100564292088
 H  -0.76584954941485162       0.47045512915600951       0.71197774484952547
 H   0.33813250754609325       -1.0217171694056011       0.43081043111330708
 H  -0.66885366253257239      -0.88565362397310610      -0.76230742051520683
 H    0.0000000000000000        0.0000000000000000        10.560402393341064
           6
           2   10.570185661315918        59.077923490250427
 C    4.9995719662588062E-002   1.9433570273590349E-002  -5.6145657920771551E-002
 H   0.54699863563179474      -0.51274152448620047       0.86368685674454160
 H   -1.3059409910751345      -0.40217994502678356       0.36928993237443070
 H    8.5167241702765351E-002 -0.56611080900603028      -0.98537152036295583
 H    7.8484625857204948E-002   1.2496400791956357       0.42091148271299739
 H    0.0000000000000000        0.0000000000000000        10.570185661315918
etc...
```


---------------------------------------------------------------------------------------
## III. OUTPUT FILES

By default, three files are created upon running PIPPy. Additional output files can be written using the options in Record 2 of the standard input.

### III.A. Standard output

The standard output contains a variety of general information, including a progress
counter for the symmetrization step. Additional debug output can be turned on by
writing to unit 1 in Recrond 2 of the standard input.

### III.B. basis.dat

The first line lists the number of atoms, atom pairs, coefficients, and terms in the expansion. Following this is a list of each term, the group it belongs to, and the exponents for each term.

### III.C. coef.dat

A list of the coefficients optimized against the training data. Each line lists the group number followed by the coefficient's value.

### III.D. Additional output

These output files are turned on using Record 2 of the standard input.

 10 : (basisinfo.dat) Write detailed basis set information, including bond distance labels\
 11 : (vtrain.dat) Compare ab intio training energies and fitted energies\
 12 : (vtest.dat) Compare ab intio test energies and fitted energies\
 20 : (xmat.dat) Write the X matrix containing the symmetrized basis\
      functions (columns) evaluated for the training data (rows)

