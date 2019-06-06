#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytfa

from cobra.io import load_matlab_model

from pytfa.io import load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from pytfa.optim.variables import DeltaG,DeltaGstd,ThermoDisplacement
from pytfa.analysis import  variability_analysis,           \
                            apply_reaction_variability,     \
                            apply_generic_variability,       \
                            apply_directionality

from cobra.io import load_matlab_model, load_json_model


from pytfa.io import import_matlab_model, load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'

# # Load the cobra_model
# cobra_model = load_matlab_model('../models/small_ecoli.mat')

# # Load reaction DB
# thermo_data = load_thermoDB('../data/thermo_data.thermodb')
# lexicon = read_lexicon('../models/small_ecoli/lexicon.csv')
# compartment_data = read_compartment_data('../models/small_ecoli/compartment_data.json')

cobra_model = load_json_model('../models/iJO1366.json')

thermo_data = load_thermoDB('../data/thermo_data.thermodb')
lexicon = read_lexicon('../models/iJO1366/lexicon.csv')
compartment_data = read_compartment_data('../models/iJO1366/compartment_data.json')

# Initialize the cobra_model
mytfa = pytfa.ThermoModel(thermo_data, cobra_model)

# Annotate the cobra_model
annotate_from_lexicon(mytfa, lexicon)
apply_compartment_data(mytfa, compartment_data)

# Initialize the cobra_model
tmodel = pytfa.ThermoModel(thermo_data, cobra_model)
tmodel.name = 'tutorial'

# Annotate the cobra_model
annotate_from_lexicon(tmodel, lexicon)
apply_compartment_data(tmodel, compartment_data)

# Set the solver
tmodel.solver = GUROBI

## TFA conversion
tmodel.prepare()
tmodel.convert(add_displacement = True)

## Info on the cobra_model
tmodel.print_info()

## Optimality
solution = tmodel.optimize()

# Apply the directionality of the solution to the cobra_model
fixed_directionality_model = apply_directionality(tmodel, solution)

# Calculate variability analysis on all continuous variables
tva_fluxes = variability_analysis(fixed_directionality_model, kind='reactions')
thermo_vars = [DeltaG,DeltaGstd,ThermoDisplacement]
tva_thermo = variability_analysis(fixed_directionality_model, kind=thermo_vars)

tight_model = apply_reaction_variability(fixed_directionality_model, tva_fluxes)
tight_model = apply_generic_variability (tight_model               , tva_thermo)

## Sample space
from pytfa.optim import strip_from_integer_variables
from pytfa.analysis import sample

continuous_model = strip_from_integer_variables(tight_model)
sampling = sample(continuous_model, 10, processes = 10)


directory = 'outputs/'
if not os.path.exists(directory):
    os.makedirs(directory)

###########################
###       Escher        ###
###########################


## Extract variable values for visualisation with Escher
from pytfa.io.viz import export_variable_for_escher

# Oppose sign to have reaction in the right direction
thermo_values = -1*sampling.median()

export_variable_for_escher(tmodel = continuous_model,
                           variable_type = ThermoDisplacement,
                           data = thermo_values,
                           filename = directory+'tutorial_sampling_median_thermo_disp.csv')


# Oppose sign to have reaction in the right direction
thermo_values = -1*sampling.mean()
export_variable_for_escher(tmodel = continuous_model,
                           variable_type = ThermoDisplacement,
                           data = thermo_values,
                           filename = directory+'tutorial_sampling_mean_thermo_disp.csv')

## Extract energy dissipation from a solution
from pytfa.analysis import calculate_dissipation
dissipation = calculate_dissipation(tmodel, solution)
dissipation.to_csv(directory+'tutorial_sampling_dissipation.csv')


###########################
###     Plotting        ###
###########################

from pytfa.io.plotting import plot_histogram
from bokeh.plotting import show, output_file

filename = directory+'tutorial_sampling_fig_{}_histogram.html'

for vid in ['LnGamma_PFK', 'LnGamma_FBA']:
    output_file(filename.format(vid))
    values = sampling[vid]
    p = plot_histogram(values)
    p.title.text = 'Samples of {}'.format(vid)
    show(p)
