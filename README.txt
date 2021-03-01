Autofit
A Fortran code for fitting permutationally invariant polynomials

Daniel R. Moberg and Ahren W. Jasper*
Argonne National Laboratory, 2021

Email: ajasper@anl.gov

---------------------------------------------------------------------------------------
I. INSTALLATION
  1. Install Fortran code

The Fortran code is located in the src/ directory, along with a compile script:
  ./compile.sh in autofit/src

Example input and output can be found in:
  autofit/examples/

  2. PIPPy

PIPPy is an optional Python wrapper for generating autofit input files.
Example input can be found in:
  autofit/examples/pippy/

Note: Mako Templates module for Python is required, see https://www.makotemplates.org/

---------------------------------------------------------------------------------------
II. INPUT FILES

Autofit requires a primary input file, here called fit.in. It also reads data from 
files containing training and test sets. These files should be located in the same 
directory the program is run and their names are specified in Record 1 of fit.in.

  1. fit.in

The Fortran input file contains a series of records, each a single line of input 
parameters.  Alternatively, input can be formatted in a more user-friendly file and 
converted with accompanying Python wrapper (located in /pippy). Below is a list of 
each record and associated parameters, with data type in parentheses.

Key: 
dp = double precision; int = integer; char * n = character*n; log = logical
[] designates number of data entries if greater than 1

Record 1: DataTrain  DataTest
training.dat  test.dat
- files (char * 10): Names of the in-sample (training) and out-of-sample (test) 
data files

Record 2: NumWrite  ITmp
nwrite   itmp[nwrite]
- nwrite (int): flag for output files to write
    = 0:  all sections are written, itmp not read
    > 0:  write output to unit 6 and each unit listed in itmp
- itmp[nwrite] (int): list of units to write to. See III. OUTPUT FILES for key.

Record 3: RangeParameter  RefEnergy
epsilon  vvref
- epsilon (dp): Range parameter for weight function
- vvref (dp): Reference energy for weight function
For more info, see Jasper and Davis, JPCA, 123(16), 3464, (2019)

Record 4: NumRanges  EnergyRanges
nc cut[1] cut[2] cut[3] ... cut[nc]
- nc (int): Number of energy ranges
- cut[nc] (dp): Values of ranges with cut[1] > cut[2] > ... > cut[nc]

Record 5: NumAtoms
natom
- natom (int): Number of atoms in system

Record 6: Symbols
symb[1] symb[2] ... symb[natom]
- symb[natom] (char * 2): Atom symbol for each atom

Record 7: AtomGroups
iagroup[1] iagroup[2] ... iagroup[natom]
- iagroup[natom] (int): Atom group for each atom

Record 8: ReadBasis, FactorOrder, TotalOrder, IMode
lreadbasis  ipow  ipowt  imode
- lreadbasis (log)
    = TRUE: read basis from basis.dat
    = FALSE: compute basis and write to basis.dat
- ipow (int): Maximum allowed order for a single factor in a term
- ipowt (int): Maximmum total order allowed for a term
- imode (int): Sets the mode for generating the PIP expansion
    = -1: Use intermolecular terms only
    =  0: Use all terms in PIP expansion
    =  1: Remove unconnected terms from basis
    =  2: Remove unconnected and intramolecular-only terms from basis
    =  3: Remove intramolecular-only terms from basis

Record 9: NumChannels
numfragchan
- numfragchan (int): Number of fragment channels to consider

Record 10: FragmentGroups
fgroup[natom]
- fgroup[natom] (int): Fragment groups in fragment channel(s)

  2. ai.all, ai.test

These files contain the training and test sets of the systems. Autofit does not 
include functionality for creating these files, which is left up to the decision 
of the user. The format for both files are the same:

Line 1: 		    nconfig
Line 2: 		    natom
Line 3:			    config#   dummy   V/cm-1
Line 4 to 4+natom: 	symb   x   y   z
Repeat lines 2-4+natom for nconfigs.

Example files can be found in each example system in:
  autofit/examples/

---------------------------------------------------------------------------------------
III. OUTPUT FILES

In general, three files are created upon running autofit. Example output files 
can be found in each example system in:
  autofit/examples/

  1. fit.out

This is the primary output file. It lists the following info on the 
expansion generation:
- Options selected (PIP size, Mode)
- System info (Atoms, Groups, Fragment Channels, Permutations up to first 100)
- Progress of symmetrization of expansion
- Info on the number of unconnected and intramolecular-only terms removed (if these options are turned on)
- Final size of basis, both terms and groups

Following the expansion generation, info on the fitting procedure is provided:
- Size of training set and number of configurations used in the fit
- Weight function parameters used
- Weighted errors found in test set for each cut specified in input
- Comparison of low energy points found in data and fit
- Out of sample test set error
- Last line contains a summary: Number of coefficients in expansion, error below cut[2], total error

  2. basis.out

The first line lists the number of atoms, atom pairs, coefficients, and terms 
in the expansion. Following this is a list of each term, the group it belongs to, 
and the exponents for each term.

  3. coef.dat

A list of the coefficients found from the fitting procedure. Each line lists the 
group number followed by the coefficient value.

  4. Units

Unit 6 (standard output): General output to fit.in
Unit 1: Extra info included in standard output
Unit 10: Detailed info on basis
Unit 11: vfit vs vai for test set
Unit 12: vfit vs vai for training set
Unit 20: Generates x matrix
