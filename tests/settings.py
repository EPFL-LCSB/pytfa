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
