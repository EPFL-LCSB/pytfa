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

class LumpGEM()
    """
    A class encapsulating the LumpGEM algorithm
    """
    def __init__(self, GEM, rcore, carbon_intake, thermo_data_path):
        """
        : param GEM: the GEM 
        : type GEM: cobra model
        : param core: list of Core reactions
        : type core: model.reactions
        : param carbon_intake: the amount of carbon atoms the cell intakes from its surrounding
        : type carbon_intake: float

        : param thermo_data_path: the path to the .thermodb database
        : type thermo_data_path : string
        """

        self._GEM = GEM

        # Extracting all reactions that lead to BBB
        self._rBBB = [rxn for rxn in GEM.reactions if "Biomass" in rxn.id]
        # Core reactions
        self._rcore = rcore
        #Core metabolites --- Is it true ? Or must mcore be user-defined ?
        self.mcore = [met for rxn in self._rcore for met in rxn.metabolites]
        # Non core reactions
        self._rncore = [rxn for rxn in _GEM.reactions if not (rxn in _rcore or rxn in _rBBB)]

        # Carbon intake
        self._C_intake = carbon_intake

        # Load reactions DB
        self._thermo_data = load_thermo_DB(thermo_data_path)


    def build_new_model():
        """
        TODO : Generate a new GEM model which will be optimized
        """
        # TODO
        self._model = Model()


    def generate_binary_variables():
        """
        Generate binary variables for each non-core reaction
        """
        self._bin_vars = {rxn : Variable(name=rxn.id, type='binary') for rxn in _rncore}

    
    def generate_constraints():
        """
        Generate carbon intake related constraints for each non-core reaction
        """
        constraints = []
        for rxn in self._rncore:
            rxn_const = Constraint(rxn.forward_variable + rxn.reverse_variable + _C_intake*_bin_vars[rxn], ub=_C_intake)
            constraints.append(rxn_const)

        self._model.add(constraints)


    def apply_thermo_constraints():
        """
        Apply thermodynamics constraints defined in thermoDB to Mcore & Rcore
        """
        # To apply the thermodynamics constraints to Rcore & Mcore only, we will remove every 
        # non core element from self._thermo_data

        # How to extract some particular metabolites / reactions from this DB ? 
        # Hard to understand its structure
        #core_thermo_data = process(self._thermo_data)

       mytfa = pytfa.ThermoModel(core_thermo_data, self._model) 

