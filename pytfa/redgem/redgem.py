#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
.. module:: redgem
   :platform: Unix, Windows
   :synopsis: RedGEM Algorithm

.. moduleauthor:: pyTFA team

Model class
"""

from pytfa.redgem.network_expansion import NetworkExpansion
from pytfa.redgem.lumpgem import LumpGEM
from .utils import remove_blocked_reactions
import yaml

class RedGEM():
    def __init__(self, gem, parameters_path, inplace=False):
        # If inplace is True, no deepcopy is performed : the modifications are applied directly onto the gem
        if inplace:
            self._gem = gem
        else:
            self._gem = gem.copy()

        # This one is used to perform the lumping
        self._source_gem = gem
        self.logger = self._gem.logger

        with open(parameters_path, 'r') as stream:
            try:
                self.params = yaml.safe_load(stream)
                print("Opened parameters file")
            except yaml.YAMLError as exc:
                print(exc)

        # If auto is activated, automatically extracts inorganics from the gem
        if "inorganic" not in self.params or self.params["inorganics"] == "auto":
            self.logger.info("Automatically computing inorganics to use")
            self.params["inorganics"] = self._extract_inorganics()

        if "force_solve" not in self.params:
            self.params["force_solve"] = False

        if "timeout" not in self.params:
            self.logger.info("Using default timeout : 3600s")
            self.params["timeout"] = 3600

        if "feasibility" not in self.params:
            self.logger.info("Using default solver feasibility : 1e-9")
            self.params["feasibility"] = 1e-9
        else:
            # numbers like 1e-9 are detected as strings by yaml module
            # to enable their use, we cast them into floats
            try:
                self.params["feasibility"] = float(self.params["feasibility"])
            except ValueError as v:
                self.logger.error(v)

        self.set_solver()

    def set_solver(self):
        if "solver" not in self.params or self.params["solver"].lower() == "auto":
            return None
        elif 'gurobi' in self.params["solver"].lower():
            solver = 'gurobi'
        elif 'cplex' in self.params["solver"].lower():
            solver = 'cplex'
        elif 'glpk' in self.params["solver"].lower():
            solver = 'glpk'
        else:
            solver = self.params["solver"]

        self._gem.solver = solver
        self._source_gem.solver = solver


    def run(self):
        # Extracting parameters
        core_subsystems = self.params["core_subsystems"]
        extracellular_system = self.params["extracellular_system"]
        biomass_rxn_ids = self.params["biomass_rxns"]

        biomass_rxns = [self._gem.reactions.get_by_id(x) for x in biomass_rxn_ids]
        main_bio_rxn = biomass_rxns[0]

        carbon_uptake = self.params["carbon_uptake"]
        growth_rate = self.params["growth_rate"]

        small_metabolites = self.params["small_metabolites"]
        cofactor_pairs = self.params["cofactor_pairs"]
        # Flatten cofactor_pairs list
        cofactors = [cofactor for pair in cofactor_pairs for cofactor in pair]
        inorganics = self.params["inorganics"]

        d = self.params["d"]
        n = self.params["n"]
        lump_method = self.params["lump_method"]

        force_solve = self.params["force_solve"]
        timeout = self.params["timeout"]
        self._gem.solver.configuration.tolerances.feasibility = self.params["feasibility"]
        self._gem.solver.configuration.tolerances.integrality = self.params["feasibility"]
        self._source_gem.solver.configuration.tolerances.feasibility = self.params["feasibility"]
        self._source_gem.solver.configuration.tolerances.integrality = self.params["feasibility"]

        self.logger.info("Computing network expansion...")
        expander = NetworkExpansion(self._gem, core_subsystems, extracellular_system,
                                    cofactors, small_metabolites, inorganics,
                                    d, n)
        reduced_gem = expander.run()
        self.logger.info("Done.")

        # Add the expansion to core reactions
        core_reactions = reduced_gem.reactions

        self.logger.info("Computing lumps...")
        lumper = LumpGEM(self._source_gem, core_reactions, self.params)
        lumps = lumper.compute_lumps(force_solve, method = lump_method)
        self.logger.info("Done.")

        self.logger.info("Create final network...")
        for rxns in lumps.values():
            reduced_gem.add_reactions(rxns)
        self.logger.info("Done.")

        reduced_gem.add_reactions(biomass_rxns)
        reduced_gem.add_reactions(lumper._exchanges)
        reduced_gem.add_reactions(lumper._transports)

        reduced_gem.objective = main_bio_rxn
        reduced_gem.reactions.get_by_id(main_bio_rxn.id).lower_bound = growth_rate

        if self.params['remove_blocked_reactions']:
            self.logger.info('Detecting blocked reactions')
            # Remove blocked reactions
            nrxn_1 = len(reduced_gem.reactions)
            self.removed_reactions = remove_blocked_reactions(reduced_gem)
            nrxn_2 = len(reduced_gem.reactions)
            self.logger.info('Removed {} blocked reaction with '
                             'FVA post-processing'.format(nrxn_1-nrxn_2))

        if main_bio_rxn.id not in reduced_gem.reactions:
            raise RuntimeError('Main Biomass reaction appears blocked')


        # For debugging purposes
        self.lumper = lumper
        main_bio_rxn.lower_bound = 0

        return reduced_gem

    def _extract_inorganics(self):
        """
        Extract inorganics from self._gem based on their formula

        :return: list of inorganics metabolites
        """

        inorganics = []
        for met in self._gem.metabolites:
            if not met.elements == {}: # Edge case
                # met is inorganic if it has 0 carbon in its formula
                if (not 'C' in met.elements) or met.elements['C'] <= 0:
                    inorganics.append(met.id)

        return inorganics



