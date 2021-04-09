"""
input file writing routines for PIPPy
"""

import os
import inpparser

def write_input_mako(job_dir_path,
                         data_train, data_test, num_write, units, #energy_units, 
                         range_parameter, ref_energy, num_ranges, energy_ranges,
                         num_atoms, symbols, atom_groups, read_basis, factor_order, total_order,
                         imode, num_channels, fragment_groups):
    """ write PIPPy input file
    """

    inp_str = inpparser.writer.pippy_input(
        data_train, data_test, num_write, units, #energy_units, 
        range_parameter, ref_energy, num_ranges, energy_ranges,
        num_atoms, symbols, atom_groups, read_basis, factor_order, total_order,
        imode, num_channels, fragment_groups)

    job_file_path = os.path.join(job_dir_path, 'input')
    with open(job_file_path, 'w') as input_file:
        input_file.write(inp_str)

def write_input(job_dir_path,
                    data_train, data_test, num_write, units, #energy_units, 
                    range_parameter, ref_energy, num_ranges, energy_ranges,
                    num_atoms, symbols, atom_groups, read_basis, factor_order, total_order,
                    imode, num_channels, fragment_groups):
    """ write PIPPy input file
    """

    inp_str = inpparser.writer.pippy_input(
        data_train, data_test, num_write, units, #energy_units, 
        range_parameter, ref_energy, num_ranges, energy_ranges,
        num_atoms, symbols, atom_groups, read_basis, factor_order, total_order,
        imode, num_channels, fragment_groups)

    job_file_path = os.path.join(job_dir_path, 'input')
    with open(job_file_path, 'w') as input_file:
        input_file.write(inp_str)
