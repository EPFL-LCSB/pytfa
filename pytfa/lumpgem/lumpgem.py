#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ..io.base import import_matlab_model, load_thermoDB
from cobra.io import load_json_model, load_yaml_model, read_sbml_model

from ..optim.variables import BinaryVariable
from ..thermo.tmodel import ThermoModel

from numpy import sum

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'


class LumpGEM:
    """
    A class encapsulating the LumpGEM algorithm
    """
    def __init__(self, path_to_model, biomass_rxns, core_subsystems, carbon_uptake, growth_rate,  thermo_data_path):
        """
        : param GEM: the GEM 
        : type GEM: cobra model

        : param biomass_rxns: list of biomass reactions
        : type biomass_rxns: [GEM.biomass_rxn.id]

        : param core_subsystems: list of Core subsystems
        : type core_subsytems: [[model.reactions]]

        : param carbon_intake: the amount of carbon atoms the cell intakes from its surrounding
        : type carbon_intake: float

        : param thermo_data_path: the path to the .thermodb database
        : type thermo_data_path : string
        """

        # Load the GEM through the appropriate cobra loading function (based on path suffix)
        self._model = self._load_model(path_to_model)

        # Extracting all reactions that lead to BBB
        self._rBBB = set([rxn for rxn in self._model.reactions if rxn.id in biomass_rxns])

        # Set containing every core reaction
        self._rcore = set([])
        # Set containing every core metabolite
        self._mcore = set([])
        for subsystem in core_subsystems:
            for rxn in subsystem:
                # Add rxn to core reactions
                self._rcore.add(rxn)
                # Add involved metabolites to core metabolites
                for met in rxn.metabolites:
                    self._mcore.add(met)

        # Non core reactions
        self._rncore = set([rxn for rxn in self._model.reactions if not (rxn in self._rcore or rxn in self._rBBB)])

        # Carbon uptake
        self._C_uptake = carbon_uptake
        # Growth rate
        self._growth_rate = growth_rate

        # TODO : solver choice
        self._solver = 'optlang-cplex'

        # Build thermo model
        self._tfa_model = self._apply_thermo_constraints(thermo_data_path, self._model)

        self._bin_vars = self._generate_binary_variables()
        self._generate_constraints()

    def _load_model(self, path_to_model):
        # Matlab
        if path_to_model[-4:] == ".mat":
            return import_matlab_model(path_to_model)

        # YAML
        if path_to_model[-4:] == ".yml":
            return load_yaml_model(path_to_model)

        # JSON
        if path_to_model[-5:] == ".json":
            return load_json_model(path_to_model)

        # SBML
        if path_to_model[-4:] == ".xml":
            return read_sbml_model(path_to_model)

    def _generate_constraints(self):
        """
        Generate carbon intake related constraints for each non-core reaction and 
        growth rate related constraints for each BBB reaction
        """
        # Carbon intake constraints
        for rxn in self._rncore:
            # rxn constrained according to the carbon uptake
            rxn_const = self._tfa_model.problem.Constraint(rxn.forward_variable +
                                                           rxn.reverse_variable +
                                                           self._C_uptake * self._bin_vars[rxn], ub=self._C_uptake)
            self._tfa_model.add_cons_vars(rxn_const)

        # Growth rate constraints
        for bio_rxn in self._rBBB:
            bio_rxn.lower_bound = self._growth_rate

    def _generate_objective(self):
        """
        Generate and add the maximization objective
        """
        # List of all binary variables
        # TODO Check the correct construction of variables
        bin_vars = [BinaryVariable(rxn.id, self._tfa_model) for rxn in self._rncore]
        # Add all binary variables to the model
        self._tfa_model.add_cons_vars(bin_vars)

        # Sum of binary variables to be maximized
        objective_sum = sum(bin_vars)
        # Set the sum as the objective function
        self._tfa_model.objective = self._tfa_model.problem.Objective(objective_sum, direction='max')

    def _apply_thermo_constraints(self, thermo_data_path, cobra_model):
        """
        Apply thermodynamics constraints defined in thermoDB to Mcore & Rcore
        """
        thermo_data = load_thermoDB(thermo_data_path)
        tfa_model = ThermoModel(thermo_data, cobra_model)
        tfa_model.name = 'Lumped Model'

        # TODO : Check what are these operations for
        # self.read_lexicon = read_lexicon()
        # compartment_data = read_compartment_data()
        # annotate_from_lexicon(tfa_model, lexicon)
        # apply_compartment_data(tfa_model, compartment_data)

        return tfa_model

    def run_optimisation(self):
        self._tfa_model.prepare()

        # Deactivate tfa computation for non-core reactions
        for ncrxn in self._rncore:
            ncrxn.thermo['computed'] = False

        self._tfa_model.convert()

        tfa_solution = self._tfa_model.optimize()
        return tfa_solution

    def lump_reaction(self, rBBB):
        """
        :param rBBB: The objective biomass reaction
        :return:
        """

        # Setting the new objective
        self._model.objective = rBBB
        # Computing TFA solution
        solution = self.run_optimisation()

        #TODO : generate lumped reaction
