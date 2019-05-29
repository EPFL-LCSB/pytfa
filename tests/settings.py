# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Initialization of globals for tests

"""
import os
import pytfa.io

from pytfa.io.enrichment import read_lexicon, annotate_from_lexicon, \
    read_compartment_data, apply_compartment_data
from cobra.test import create_test_model

############
# SETTINGS #
############

# Objective value of the Matlab solution
objective_value = 0.810997250260066

######## End of Settings ########


this_directory = os.path.dirname(os.path.realpath(__file__))

# Load the thermo database
thermo_data = pytfa.io.load_thermoDB \
                                    (this_directory \
                                     + '/../data/thermo_data.thermodb')

# Load the cobra_model
cobra_model = pytfa.io.import_matlab_model \
                                    (this_directory \
                                     + '/../models/small_ecoli.mat')

# Make your computations on it
tmodel = pytfa.ThermoModel(thermo_data, cobra_model)
tmodel.solver = 'optlang-glpk'
tmodel.prepare()
tmodel.convert()


# Small model for simpler tests

small_model = create_test_model('textbook')

# Make your computations on it
# tmodel = pytfa.ThermoModel(thermo_data, cobra_model)

lexicon = read_lexicon(this_directory+'/../models/iJO1366/lexicon.csv')
compartment_data = read_compartment_data(this_directory+'/../models/iJO1366/compartment_data.json')

# Initialize the cobra_model
small_tmodel = pytfa.ThermoModel(thermo_data, small_model)

# Annotate the cobra_model
annotate_from_lexicon(small_tmodel, lexicon)
apply_compartment_data(small_tmodel, compartment_data)

small_tmodel.solver = 'optlang-glpk'
small_tmodel.prepare()
small_tmodel.convert()