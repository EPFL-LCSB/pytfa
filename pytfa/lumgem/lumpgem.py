#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. module:: lumpgem
   :platform: Unix, Windows
   :synopsis: LumpGEM Algorithm

.. moduleauthor:: pyTFA team

Model class
"""

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

# Temporary class parameter, until the load-from-paramfile feature is added
class LumpGEM()
    def __init__(self, GEM, core, carbon_intake):
        """
        : type GEM cobra model
        : param GEM the GEM 
        : type core model.reactions
        : param core list of Core reactions
        """

        self._GEM = GEM

        # Extracting all reactions that lead to BBB
        self._rBBB = [rxn for rxn in GEM.reactions if "Biomass" in rxn.id]
        # Core reactions
        self._rcore = core
        # Non core reactions
        self._rncore = [rxn for rxn in _GEM.reactions if not (rxn in _rcore or rxn in _rBBB)]

        # Carbon intake
        self._C_intake = carbon_intake


    def build_new_model():
        return

    def generate_binary_variables():
        self._bin_vars = {rxn : Variable(name=rxn.id, type='binary') for rxn in _rncore}
    
    def generate_constraints():
        constraints = []
        for rxn in model.reactions:
            rxn_const = Constraint(rxn.forward_variable + rxn.reverse_variable + _C_intake*_bin_vars[rxn], ub=_C_intake)
            constraints.append(rxn_const)


        _model.add(constraints)
