#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytfa

from pytfa.io import import_matlab_model, load_thermoDB

from pytfa.optim.variables import DeltaG,DeltaGstd,ThermoDisplacement
from pytfa.analysis import  variability_analysis,           \
                            apply_reaction_variability,     \
                            apply_generic_variability,       \
                            apply_directionality

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'

# Load the model
cobra_model = import_matlab_model('../models/glycolysis.mat')

# Load reaction DB
thermo_data = load_thermoDB('../data/thermo_data.thermodb')

# Initialize the model
tmodel = pytfa.ThermoModel(thermo_data, cobra_model)
tmodel.normalize_reactions()

# Set the solver
tmodel.solver = CPLEX

## TFA conversion
tmodel.prepare()
tmodel.convert(add_displacement = True)

## Info on the model
tmodel.print_info()

## Optimality
solution = tmodel.optimize()

# Apply the directionality of the solution to the model
fixed_directionality_model = apply_directionality(tmodel, solution)

# Calculate variability analysis on all continuous variables
tva_fluxes = variability_analysis(fixed_directionality_model, kind='reactions')
tight_model = apply_reaction_variability(fixed_directionality_model, tva_fluxes)


thermo_vars = [DeltaG,DeltaGstd,ThermoDisplacement]
tva_thermo = variability_analysis(tight_model, kind=thermo_vars)

tight_model = apply_generic_variability (tight_model, tva_thermo)

## Sample space
from pytfa.optim import strip_from_integer_variables
from pytfa.analysis import sample

continuous_model = strip_from_integer_variables(tight_model)
sampling = sample(continuous_model, 10000, processes = 4)

directory = 'outputs/'
if not os.path.exists(directory):
    os.makedirs(directory)

###########################
###       Escher        ###
###########################


## Extract variable values for visualisation with Escher
from pytfa.io.viz import export_variable_for_escher, export_reactions_for_escher

# Oppose sign to have reaction in the right direction
thermo_values = -1*sampling.median()

export_variable_for_escher(tmodel = continuous_model,
                           variable_type = ThermoDisplacement,
                           data = thermo_values,
                           filename = directory+'glycolysis_median_thermo_disp.csv')

# Oppose sign to have reaction in the right direction
thermo_values = -1*sampling.mean()
export_variable_for_escher(tmodel = continuous_model,
                           variable_type = ThermoDisplacement,
                           data = thermo_values,
                           filename = directory+'glycolysis_mean_thermo_disp.csv')

export_reactions_for_escher(tmodel,
                           sampling,
                           directory+'glycolysis_avg_fluxes.csv')


###########################
###     Plotting        ###
###########################

from pytfa.io.plotting import plot_histogram
from bokeh.plotting import show, output_file

filename = directory+'glycolysis_fig_{}_histogram.html'

for rid in ['PGK', 'FBP', 'GLCptspp', 'GAPD' ]:
    vid = 'LnGamma_' + rid
    output_file(filename.format(vid))
    values = sampling[vid]
    p = plot_histogram(values)
    p.title.text = 'Log Thermodynamic displacement of {}'.format(rid)
    show(p)