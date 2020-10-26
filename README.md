# autofit
## automated PIP fitting code

## Authors
Ahren W. Jasper\
Daniel R. Moberg

## Functionality

## Contact
Ahren Jasper [ajasper@anl.gov]

## Install
./compile.sh in autofit/src

For pyfit wrapper input generator:\
   Mako Templates module for Python is required\
   see [https://www.makotemplates.org/]

## Fortran input file
The Fortran input file contains a series of records, each a single line of input parameters (see example fit.in files)
Alternatively, input can be formatted in more user-friendly file and converted with accompanying Python wrapper (located in /pyfit)
Below are listed each record and associated parameters, with data type in parantheses

**Record 1**: TrainingFile, TestFile\
training.dat  test.dat
- files (character * 10) Names of the training and test data files

**Record 2**: RangeParam and RefEnergy\
epsilon enref
- epsilon (dp) Range parameter for weight function
- enref (dp) Reference energy for weight function (see Jasper and Davis, JPCA, 123(16), 3464, (2019))

**Record 3**: EnergyRanges\
nc cut(1) cut(2) cut(3) ... cut(nc)
- nc (int) Number of energy ranges
- cut() (dp) Values of ranges with cut(1) > cut(2) > ... > cut(nc)

**Record 4**: NumAtoms\
natom (int)
- Number of atoms in system

**Record 5**: Symbols\
symb(natom) (char * 2)
- Atom symbol for each atom

**Record 6**: AtomGroups\
iagroup(natom) (int)
- Atom group for each atom

**Record 7**: FactorOrder, TotalOrder, ReadBasis, ExpType, NumFragChannels\
ipow ipowt lreadbasis exptype numfragchan
- ipow (int) Maximum allowed order for a single factor in a term
- ipowt (int) Maximmum total order allowed for a term
- lreadbasis (logical)
    = TRUE, read basis from basis.dat
    = FALSE, compute basis and write to basis.dat
- exptype (int)
    - 0: Use all terms in PIP expansion
    - 1: Remove unconnected terms and intramolecular terms from basis
    - 2: Remove unconnected terms
    - 3: Remove intra-only terms
    - -1: Use only intermolecular terms
- numfragchan (int) Number of possible fragment channels to consider for Record 8

**Record 8**: Fragment Groups\
fgroup(natom) (int)
- Molecular groups in the product system(s).

