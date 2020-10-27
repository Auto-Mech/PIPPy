"""
Writes input for PIPPy
"""

import os
from mako.template import Template


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def pippy_input(data_sets, energy_units, range_parameter, ref_energy, 
        num_ranges, energy_ranges, 
        num_atoms, symbols, atom_groups, factor_order, total_order,
        read_basis, exp_type, num_channels, fragment_groups):
    """ writes the PIPPy input file for each instance
    """

    # Set the dictionary for the PIPPy input file
    fill_vals = {
        "DataSets": data_sets,
        "EnergyUnits": energy_units,
	"RangeParameter": range_parameter,
	"RefEnergy": ref_energy,
	"NumRanges": num_ranges,
        "EnergyRanges": energy_ranges,
        "NumAtoms": num_atoms,
        "Symbols": symbols,
        "AtomGroups": atom_groups,
        "FactorOrder": factor_order,
        "TotalOrder": total_order,
        "ReadBasis": read_basis,
	"ExpansionType": exp_type,
	"NumChannels": num_channels,
        "FragmentGroups": fragment_groups
    }

    # Set template name and path for the PIPPy input file
    template_file_name = 'pippy.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)
    print("Template path",template_file_path)

    # Build the 1dmin input string
    input_str = Template(filename=template_file_path).render(**fill_vals)

    return input_str
