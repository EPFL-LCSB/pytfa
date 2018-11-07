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
    def __init__(self, GEM, core):
        """
        : type GEM cobra model
        : param GEM the GEM 
        : type core model.reactions
        : param core list of Core reactions
        """

        self._GEM = GEM
        self._rcore = core
        #self._mcore = [metab for metab in self._rcore.metabolites.keys()]

        # Extracting all reactions that lead to BBB
        self.BBBreactions = [rxn for rxn in GEM.reactions if "Biomass" in rxn.id]

