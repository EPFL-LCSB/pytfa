"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Some tools used by pyTFA


"""


import re


# Regular Expression to get the atoms in a reaction. Compile it once for all


def camel2underscores(name):
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def varnames2ids(tmodel, variables):
    return [tmodel._var_dict[x].id for x in variables]