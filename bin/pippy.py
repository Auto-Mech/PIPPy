""" Sets up PIPPy
"""

import os
import sys
import subprocess
import ioformat
import pippy_io


# Set paths
DRIVE_PATH = os.getcwd()

# Parse the input file into a string
INP_STR = ioformat.pathtools.read_file(
    DRIVE_PATH, 'pippy.inp', remove_comments='#', remove_whitespace=True)

print('Parsing initial input...')
TRAIN_DCT, FFORM_DCT, EXEC_DCT = pippy_io.parser.input_file(INP_STR)

# Write input file
print('Writing input file to input ...')
PIP_INP_STR = pippy_io.writer.input_file(
    data_train=TRAIN_DCT['DataTrain'],
    data_test=TRAIN_DCT['DataTest'],
    num_write=TRAIN_DCT['NumWrite'],
    units=TRAIN_DCT['Units'],
    range_parameter=TRAIN_DCT['RangeParameter'],
    ref_energy=TRAIN_DCT['RefEnergy'],
    num_ranges=TRAIN_DCT['NumRanges'],
    energy_ranges=TRAIN_DCT['EnergyRanges'],
    num_atoms=FFORM_DCT['NumAtoms'],
    symbols=FFORM_DCT['Symbols'],
    atom_groups=FFORM_DCT['AtomGroups'],
    read_basis=FFORM_DCT['ReadBasis'],
    factor_order=FFORM_DCT['FactorOrder'],
    total_order=FFORM_DCT['TotalOrder'],
    imode=FFORM_DCT['IMode'],
    num_channels=FFORM_DCT['NumChannels'],
    fragment_groups=FFORM_DCT['FragmentGroups']
)
ioformat.pathtools.write_file(PIP_INP_STR, DRIVE_PATH, 'input')

# Run Fortran code with new input file if desired
if EXEC_DCT['UseCL'] == 'T':
    print('Running Fortran code ...')
    print(EXEC_DCT['CommandLine'])
    RESULT = subprocess.call(EXEC_DCT['CommandLine'], shell=True)
    if RESULT == 0:
        print('Finished')
    else:
        print('Error: Command Failed')
        sys.exit()
