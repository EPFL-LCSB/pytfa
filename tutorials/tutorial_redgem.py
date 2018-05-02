#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytfa

from optlang.exceptions import SolverError
from cobra.core.model import SolverNotFound
from cobra.flux_analysis import flux_variability_analysis

from pytfa.io import import_matlab_model, load_thermoDB

cobra_model = import_matlab_model('../models/full_ecoli.mat')

# Load reaction DB
thermo_data = load_thermoDB('../data/thermo_data.thermodb')

mytfa = pytfa.ThermoModel(thermo_data, cobra_model)
mytfa.solver = 'optlang-glpk'

## FBA
fba_solution = mytfa.optimize()
fba_value = fba_solution.f
fva = flux_variability_analysis(mytfa)

## TFA conversion
mytfa.prepare()
mytfa.convert()

## Info on the cobra_model
mytfa.print_info()

## Optimality
tfa_solution = mytfa.optimize()


#------------------ Start WIP ------------------------#

from pytfa.reduction import reduce

reduction_parameters = load_params('reduction_params.py')

red_tfa = reduce(mytfa,reduction_parameters)

red_tfa.print_info()

red_tfa.optimize()