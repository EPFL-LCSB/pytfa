#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. module:: redgem
   :platform: Unix, Windows
   :synopsis: RedGEM Algorithm

.. moduleauthor:: pyTFA team

Model class
"""

import network_expansion
import lumpgem
import yaml

class RedGEM():
    def __init__(self, gem, inplace, parameters_path):
        # If inplace is True, no deepcopy is performed : the modifications are applied directly onto the gem
        if inplace:
            self._gem = gem
        else:
            self._gem = gem.copy()

        with open(parameters_path, 'r') as stream:
            try:
                self.params = yaml.safe_load(stream))
                print("Opened parameters file")
            except yaml.YAMLError as exc:
                print(exc)

    def run():
        # Extracting parameters
        core_subsystems = self.params["core_subsystems"]
        subsystem_names = self.params["subsystem_names"]
        extracellular_system = self.params["extracellular_system"]
        biomass_rxns = self.params["biomass_rxns"]

        carbon_uptake = self.params["carbon_uptake"]
        growth_rate = self.params["growth_rate"]

        small_metabolites = self.params["small_metabolites"]
        inorganics = self.params["inorganics"]

        d = self.params["d"]
        n = self.params["n"]

        timeout = self.params["timeout"]

        print("Computing network expansion...")
        expander = NetworkExpansion(self._gem, core_subsystems, carbon_uptake, cofactor_pairs, small_metabolites, inorganics, d, extracellular_system, subsystem_names, n)
        reduced_gem = expander.run()
        print("Done.")

        print("Computing lumps...")
        lumper = LumpGEM(reduced_gem, biomass_rxns, core_subsystems, carbon_uptake, growth_rate, timeout)
        lumps = lumper.run()
        print("Done.")
        return lumps

