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
import yaml

class RedGEM():
    def __init__(self, gem, parameters_path, inplace=False):
        # If inplace is True, no deepcopy is performed : the modifications are applied directly onto the gem
        if inplace:
            self._gem = gem
        else:
            self._gem = gem.copy()

        with open(parameters_path, 'r') as stream:
            try:
                self.params = yaml.safe_load(stream)
                print("Opened parameters file")
            except yaml.YAMLError as exc:
                print(exc)

        # If auto is activated, automatically extracts inorganics from the gem
        if self.params["inorganics"] == "auto":
            self.params["inorganics"] = self._extract_inorganics()

    def run(self):
        # Extracting parameters
        core_subsystems = self.params["core_subsystems"]
        extracellular_system = self.params["extracellular_system"]
        biomass_rxns = self.params["biomass_rxns"]

        carbon_uptake = self.params["carbon_uptake"]
        growth_rate = self.params["growth_rate"]

        small_metabolites = self.params["small_metabolites"]
        cofactor_pairs = self.params["cofactor_pairs"]
        inorganics = self.params["inorganics"]

        d = self.params["d"]
        n = self.params["n"]

        timeout = self.params["timeout"]

        print("Computing network expansion...")
        expander = NetworkExpansion(self._gem, core_subsystems, extracellular_system,
                                    cofactor_pairs, small_metabolites, inorganics,
                                    d, n)
        reduced_gem = expander.run()
        print("Done.")

        print("Computing lumps...")
        lumper = LumpGEM(reduced_gem, biomass_rxns, core_subsystems, carbon_uptake, growth_rate, timeout)
        lumps = lumper.compute_lumps()
        print("Done.")
        return lumps

    def _extract_inorganics(self):
        """
        Extract inorganics from self._gem based on their formula

        :return: list of inorganics metabolites
        """

        inorganics = []
        for met in self._gem.metabolites:
            if not met.elements == {}: # Edge case
                # met is inorganic if it has 0 carbon in its formula
                if 'C' in met.elements and met.elements['C'] > 0:
                    inorganics.append(met)

        return inorganics



