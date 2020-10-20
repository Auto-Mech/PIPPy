# autofit
## automated PIP fitting code

## Authors
Ahren W. Jasper
Daniel R. Moberg

## Functionality

## Contact
Ahren Jasper [ajasper@anl.gov]

# Install
./compile.sh in autofit/src

For pyfit wrapper input generator:\
   Mako Templates module for Python is required\
   see [https://www.makotemplates.org/]

# Input Files
## Fortran input file
The Fortran input file contains a series of records, each a single line of input parameters (see example fit.in files)
Alternatively, input can be formatted in more user-friendly file and converted with accompanying Python wrapper (located in /pyfit)
Below are listed each record and associated parameters, with data type in parantheses

**Record 1**: TrainingFile, TestFile\
training.dat  test.dat (character * 10)\
Names of the training and test data files

**Record 2**: RangeParam and RefEnergy\
epsilon enref (dp)\
-Range parameter and reference energy for weight function (see Jasper and Davis, JPCA, 123(16), 3464, (2019))

**Record 3**: EnergyRanges\
nc cut(1) cut(2) cut(3) ... cut(nc) (dp)\
-Number of energy ranges (nc) and their values with cut(1) > cut(2) > ... > cut(nc)

**Record 4**: NumAtoms\
natom (int)\
-Number of atoms in system

**Record 5**: Symbols\
symb(natom) (char * 2)\
-Atom symbol for each atom

**Record 6**: AtomGroups\
iagroup(natom) (int)\
-Atom group for each atom

**Record 7**: FactorOrder, TotalOrder, ReadBasis, ExpType, NumFragChannels\
ipow ipowt lreadbasis exptype numfragchan (int, int, log, int, int)\
-Maximum allowed order for each factor and total order of the term\
-lreadbasis (logical)\
    = TRUE, read basis from basis.dat\
    = FALSE, compute basis and write to basis.dat\
-exptype (logical)\
    = TRUE, determine disconnected atom groups to remove from basis\
    = FALSE, do not consider disconnected groups

**Record 8**: Fragment Groups\
fgroup(natom) (int)\
-Molecular groups in the product system(s).

## ai.all, ai.test

