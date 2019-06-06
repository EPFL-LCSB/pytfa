#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import errno
import pytfa

from pytfa.io import import_matlab_model, load_thermoDB

from pytfa.optim.variables import DeltaG,DeltaGstd,ThermoDisplacement
from pytfa.analysis import  variability_analysis,           \
                            apply_reaction_variability,     \
                            apply_generic_variability,       \
                            apply_directionality

from cobra.flux_analysis.variability import flux_variability_analysis

from math import log

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'

# Load the cobra_model
cobra_model = import_matlab_model('../models/small_ecoli.mat')

# Load reaction DB
thermo_data = load_thermoDB('../data/thermo_data.thermodb')

# Initialize the cobra_model
tmodel = pytfa.ThermoModel(thermo_data, cobra_model)

# Set the solver
tmodel.solver = GLPK

## TFA conversion
tmodel.prepare()
tmodel.convert(add_displacement = True)

## Info on the cobra_model
tmodel.print_info()

## Optimality
solution = tmodel.optimize()

# Calculate variability analysis on all continuous variables
fva_fluxes = flux_variability_analysis(cobra_model)
tva_fluxes = variability_analysis(tmodel, kind='reactions')

# Add more specific concentration data
def apply_concentration_bound(met, lb, ub):
    the_conc_var = tmodel.log_concentration.get_by_id(met)
    # Do not forget the variables in the model are logs !
    the_conc_var.ub = log(ub)
    the_conc_var.lb = log(lb)

apply_concentration_bound('atp_c', lb=1e-3, ub=1e-2)
apply_concentration_bound('adp_c', lb=4e-4, ub=7e-4)
apply_concentration_bound('amp_c', lb=2e-4, ub=3e-4)

tmodel.optimize()
# Perform variability analysis again
tva_fluxes_lc = variability_analysis(tmodel, kind='reactions')

###########################
###     Plotting        ###
###########################

from pytfa.io.plotting import plot_fva_tva_comparison
from bokeh.plotting import show, output_file
from bokeh.layouts import column

try:
    os.mkdir('outputs')
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass

output_file('outputs/va_comparison.html')
p1 = plot_fva_tva_comparison(fva_fluxes, tva_fluxes)
p2 = plot_fva_tva_comparison(fva_fluxes, tva_fluxes_lc)
c = column(p1,p2)
show(c)