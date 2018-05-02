# -*- coding: utf-8 -*-
"""
.. module:: pyTFA
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Compare two LP files, representing constraints and variables of a cobra_model

Usage: python lpdiff.py file1.py file2.py


"""


import sys
import re

FLOAT_PRECISION = 10**-3

# Regular expressions database, used to parse the data from LP files
regex = {
    'name':re.compile(r'\\Problem name: (.+)'),
    'obj' :re.compile(r'obj: (.+)'),
    'cons':re.compile(r'([0-9A-Za-z\-_]+): (.+)([=<>]+) ([\-e0-9\\.]+)'),
    'expr':re.compile(r' *([\-+]? *[0-9\\.]*(?:e-?[0-9]+)?[ \*]+|-|)([A-Za-z0-9][A-Za-z0-9_\-]+) +'),
    'number':re.compile(r'[0-9]'),
    'bounds':re.compile(r'([0-9e\.\-\+]+) +<= +(.+) +<= +([0-9e\.\-\+]+)'),
    'bin':re.compile(r'[A-Za-z0-9\-_]+'),
    'reverse':re.compile(r'(.+_reverse)_[a-f0-9]{5}')
}


def var_cleanup(name):
    """ Standardize a variable name, to make comparison between models easier

    :params str name: The variable name to standardize

    :returns: The standardized name
    :rtype: str


    """
    global regex
    # Remove spaces around the name (at the end or the beginning)
    name = name.strip()
    # Reverse flux
    if name[:2] == 'R_':
        name += '_reverse'
        name = name[2:]
    # 'R_' or 'F_'
    if name[:2] in ['F_', 'M_']:
        name = name[2:]
    # Reverse fluxes by COBRApy
    test = regex['reverse'].match(name)
    if test:
        name = test.groups()[0]
    return name


def float_cleanup(val):
    """ Convert an extracted string to a float

    :param str val: The string to convert to a float

    :returns: The float corresponding to the value given
    :rtype: float


    """
    global regex
    val = val.replace(' ', '').replace('*', '')
    if not regex['number'].search(val):
        val += '1'
    return float(val)


def process_expr(expr):
    """ Parse the string of an expression

    :param str expr: The string of the expression to parse

    :returns: A dictionnary with the variables' names as keys and their
    coefficients as values

    :rtype: dict(float)

    """
    global regex
    items = regex['expr'].findall(expr)
    res = {}
    for item in items:
        met = var_cleanup(item[1])
        if met in res:
            print(met)
            raise Exception('met already in res ! ')
        res[met] = float_cleanup(item[0])
    return res


def parse_file(path):
    """ Parse the content of a file and get the corresponding cobra_model

    :param path:
    :returns: A dictionnary representing the cobra_model contained in the file
    :rtype: dict

    """
    global regex

    with open(path) as file:
        model = dict()
        mode = 'init'
        for line in file.readlines():
            # Remove beginning and ending spaces
            line = line.strip()

            if line == '': # Skip empty lines
                continue

            if mode == 'init':
                test = regex['name'].match(line)
                if test:
                    model['name'] = test.groups()[0]
                    continue
                if line == 'Maximize':
                    mode = 'obj'
                    continue
            elif mode == 'obj':
                test = regex['obj'].match(line)
                if test:
                    model['obj'] = var_cleanup(test.groups()[0])
                    # FIXME: Multiple objectives and coefficients
                    continue
                if line == 'Subject To':
                    mode = 'cons'
                    model['cons'] = {}
                    continue
            elif mode == 'cons':
                test = regex['cons'].match(line)
                if test:
                    cons_name = var_cleanup(test.groups()[0])
                    cons_expr = test.groups()[1]
                    cons_type = test.groups()[2]
                    cons_eq   = test.groups()[3]
                    cons_eq = float(cons_eq.replace(' ',''))
                    cons_expr = process_expr(cons_expr)
                    if cons_name in model['cons']:
                        print(model['cons'][cons_name], cons_expr)
                        print(cons_name)
                        raise Exception('Duplicate constraint name !')
                    model['cons'][cons_name]  = dict()
                    model['cons'][cons_name]['expr'] = cons_expr

                    model['cons'][cons_name]['_lb'] = cons_eq
                    model['cons'][cons_name]['_ub'] = cons_eq
                    if cons_type[0] == '<':
                        model['cons'][cons_name]['_lb'] = float('-inf')
                    if cons_type[0] == '>':
                        model['cons'][cons_name]['_ub'] = float('+inf')
                    continue
                #TODO : More complex constraints (both sides) ?
                if line == 'Bounds':
                    mode = 'bounds'
                    model['bounds'] = {}
                    continue
            elif mode == 'bounds':
                test = regex['bounds'].match(line)
                if test:
                    bound_name = var_cleanup(test.groups()[1])
                    bound_ub = float_cleanup(test.groups()[2])
                    bound_lb = float_cleanup(test.groups()[0])
                    if bound_name in model['bounds']:
                        raise Exception('Duplicate Bound Name !')
                    model['bounds'][bound_name] = {
                        'lb': bound_lb,
                        'ub': bound_ub
                    }
                    continue
                if line == 'Binaries':
                    mode = 'bin'
                    model['bin'] = {}
                    continue
            elif mode == 'bin':
                bins = regex['bin'].findall(line)
                for the_bin in bins:
                    model['bin'][var_cleanup(the_bin)] = True
                if len(bins) > 0:
                    continue
                if line == 'End':
                    continue
            print(mode)
            raise Exception('Unrecognized line : ' + line)


    return model


# Now, compare them !
def compare(model, stack = None, n = 0):
    """ Deep comparison of two dictionnaries

    :param list(dict) model: A list of two dictionnaries to compare
    :param str stack: `Optional` The path where we are currently (a
        concatenation of the keys of the dictionnaries). If not set, the
        function assumes we're at the root of the dictionnaries
    :param int n: `Optional` The number of differences encountered so far. If
        not set, 0 is assumed.

    :returns: The number of differences found between the two dictionnaries
    :rtype: int

    This function will call itself recursively every time it encounters a
    dictionnary compare the content of each dictionnary. If a difference is
    found, some information is printed to STDOUT


    """
    # Make sure all keys in the first cobra_model are equal to those in the second
    for item in model[0]:
        if item not in model[1]:
            print((stack + '.' if stack else '') + item + ' in 1 but not in 2')
            n += 1
            # from IPython.thermo.debugger import Tracer
            # Tracer()()
        elif isinstance(model[0][item],dict):
            n += compare([model[0][item], model[1][item]],
                          (stack  + '.' if stack else '') + item)
        elif isinstance(model[0][item],float):
            if abs(model[0][item] - model[1][item]) > FLOAT_PRECISION:
                print((stack  + '.' if stack else '') + item + ' different :')
                print(model[0][item], ' / ', model[1][item])
                n += 1
        else:
            if model[0][item] != model[1][item]:
                print((stack + '.' if stack else '') + item + ' different :')
                print(model[0][item], ' / ', model[1][item])
                n += 1
    # Are there keys in the second cobra_model not in the first ?
    for item in model[1]:
        if item not in model[0]:
            print((stack + '.' if stack else '') + item + ' in 2 but not in 1')
            n += 1
    return n

if __name__ == '__main__':
    # Print usage information if needed
    if len(sys.argv) < 2:
        print('Usage : python ' + __file__  + ' file1.lp file2.lp')
        sys.exit(1)

    mytfa = [None] * 2
    for i in range(2):
        mytfa[i] = parse_file(sys.argv[i + 1])

    # Now, let's compare !
    print(compare(mytfa, ''))
