""" Parses the PIPPy input file
"""

import sys
import ioformat


# Required and Supported keywords for the input
SUPPORTED_KEYS = {
    'training_data': ('DataTrain', 'DataTest', 'NumWrite', 'Units',
                      'RangeParameter', 'RefEnergy',
                      'NumRanges', 'EnergyRanges'),
    'functional_form': ('NumAtoms', 'Symbols', 'AtomGroups',
                        'ReadBasis', 'FactorOrder', 'TotalOrder',
                        'IMode', 'NumChannels', 'FragmentGroups'),
    'fortran_execution': ('UseCL', 'CommandLine')
}
REQUIRED_KEYS = {
    'training_data': ('DataTrain', 'DataTest', 'NumWrite',
                      'RangeParameter', 'RefEnergy',
                      'NumRanges', 'EnergyRanges'),
    'functional_form': ('NumAtoms', 'Symbols', 'AtomGroups',
                        'ReadBasis', 'FactorOrder', 'TotalOrder',
                        'IMode', 'NumChannels', 'FragmentGroups'),
    'fortran_execution': ('UseCL',)
}

# Defaults for keys that may be missing
DEFAULT_DCT = {
    'Units': ''
}


# Parse out the input string
def input_file(inp_str):
    """ Parse the input string
    """

    # Parse the sections of the input into keyword-val dictionaries
    train_block = ioformat.ptt.symb_block(inp_str, '$', 'training_data')
    fform_block = ioformat.ptt.symb_block(inp_str, '$', 'functional_form')
    exec_block = ioformat.ptt.symb_block(inp_str, '$', 'fortran_execution')

    train_dct = ioformat.ptt.keyword_dct_from_block(
        train_block[1], formatvals=False)
    fform_dct = ioformat.ptt.keyword_dct_from_block(
        fform_block[1], formatvals=False)
    exec_dct = ioformat.ptt.keyword_dct_from_block(
        exec_block[1], formatvals=False)

    # Set defaults (maybe use fancy version later if more defaults can be set)
    if 'Units' not in train_dct:
        train_dct['Units'] = DEFAULT_DCT['Units']

    # Check that the dictionaries are built correctly
    _check_dcts(train_dct, fform_dct, exec_dct)

    return train_dct, fform_dct, exec_dct


def _check_dcts(train_dct, fform_dct, exec_dct):
    """ Assess if the dicts are build correctly
    """

    chk_info = zip((train_dct, fform_dct, exec_dct),
                   ('training_data', 'functional_form', 'fortran_execution'))
    for dct, name in chk_info:
        # Assess if a required section was defined in the input
        if dct is None:
            print('Section: {} not defined in input'.format(name))
            sys.exit()
        else:
            # Get the keywords that the user defined in the input
            defined_keys = set(dct.keys())

            # Check if only supported keys defined
            supp_keys = set(SUPPORTED_KEYS[name])
            unsupported_defined_keys = defined_keys - supp_keys
            if unsupported_defined_keys:
                print('Unsupported keywords given in section {}'.format(name))
                for key in unsupported_defined_keys:
                    print(key)
                sys.exit()

            # Check if all required keys defined
            req_keys = set(REQUIRED_KEYS[name])
            undefined_required_keys = req_keys - defined_keys
            if undefined_required_keys:
                print('Required keywords not given in section {}'.format(name))
                for key in undefined_required_keys:
                    print(key)
                sys.exit()

    # Need a secondary check on the exec dct
    if exec_dct['UseCL'] == 'T':
        if 'CommandLine' not in exec_dct:
            print('CommandLine must be given in fortran_execution ',
                  'if UseCL set to T')
            sys.exit()
