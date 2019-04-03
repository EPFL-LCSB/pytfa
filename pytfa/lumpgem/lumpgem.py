#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ..io.base import import_matlab_model, load_thermoDB
from cobra.io import load_json_model, load_yaml_model, read_sbml_model
from cobra import Reaction
from ..optim.utils import symbol_sum

from pytfa.optim.variables import ReactionVariable, BinaryVariable, get_binary_type
from pytfa.optim.constraints import ReactionConstraint
from ..thermo.tmodel import ThermoModel

from pytfa.io import read_compartment_data, apply_compartment_data, \
                     read_lexicon, annotate_from_lexicon

from numpy import sum

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'


class MyVariableClass(ReactionVariable, BinaryVariable):
    prefix = 'VC_'

    def __init__(self, reaction, **kwargs):
        ReactionVariable.__init__(self, reaction,
                                  type=get_binary_type(),
                                  **kwargs)


# Define a new constraint type:
class MyConstraintClass(ReactionConstraint):
    prefix = 'CC_'


class LumpGEM:
    """
    A class encapsulating the LumpGEM algorithm
    """
    def __init__(self, path_to_model, biomass_rxns, core_subsystems, carbon_uptake, growth_rate,  thermo_data_path, auxiliary_data):
        """
        : param GEM: the GEM 
        : type GEM: cobra model

        : param biomass_rxns: list of biomass reactions
        : type biomass_rxns: [GEM.biomass_rxn.id]

        : param core_subsystems: list of Core subsystems names
        : type core_subsystems: [string]

        : param carbon_intake: the amount of carbon atoms the cell intakes from its surrounding
        : type carbon_intake: float

        : param growth_rate: theoretical maximum specific growth rate in 1/hr units
        : type growth_rate: float

        : param thermo_data_path: the path to the .thermodb database
        : type thermo_data_path : string
        """

        # Load the GEM through the appropriate cobra loading function (based on path suffix)
        model = self._load_model(path_to_model)
        # Build thermo model
        self._tfa_model = self._apply_thermo_constraints(thermo_data_path, model, auxiliary_data)

        # Set containing every BBB reaction
        self._rBBB = []
        # Set containing every core reaction
        self._rcore = []
        # Set containing every core metabolite
        self._mcore = []
        # Set containing every non-core reaction
        self._rncore = []

        # For each reaction
        for rxn in self._tfa_model.reactions:
            # If it's a BBB reaction
            # TODO : make it possible to use keywords to match BBB reactions, rather than IDs
            if rxn.id in biomass_rxns:
                self._rBBB.append(rxn)
            # If it's a core reaction
            elif rxn.subsystem in core_subsystems:
                self._rcore.append(rxn)
                for met in rxn.metabolites:
                    self._mcore.append(met)
            # If it's neither BBB nor core, then it's non-core
            else:
                self._rncore.append(rxn)

        # Carbon uptake
        self._C_uptake = carbon_uptake
        # Growth rate
        self._growth_rate = growth_rate

        # TODO : solver choice
        self._solver = 'optlang-cplex'

        # lumpgen binary variables to deactivate non-core reactions. The reaction is deactivated when the value of
        # the variable is 1
        self._activation_vars = {rxn: self._tfa_model.add_variable(kind=MyVariableClass,
                                                                   hook=rxn,
                                                                   lb=0,
                                                                   ub=1,
                                                                   queue=False)
                                 for rxn in self._rncore}

        self._generate_carbon_constraints()
        self._generate_objective()

    def _load_model(self, path_to_model):
        # MATLAB
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

    def _apply_thermo_constraints(self, thermo_data_path, cobra_model, auxiliary_data_path):
        """
        Apply thermodynamics constraints defined in thermoDB to Mcore & Rcore
        """
        thermo_data = load_thermoDB(thermo_data_path)
        tfa_model = ThermoModel(thermo_data, cobra_model)
        tfa_model.name = 'Lumped Model'

        # TODO : Improve management of auxiliary data paths
        if auxiliary_data_path[-1] != '/':
            auxiliary_data_path += '/'
        lexicon = read_lexicon(auxiliary_data_path+'lexicon.csv')
        compartment_data = read_compartment_data(auxiliary_data_path+'compartment_data.json')

        annotate_from_lexicon(tfa_model, lexicon)
        apply_compartment_data(tfa_model, compartment_data)

        return tfa_model

    def _generate_carbon_constraints(self):
        """
        Generate carbon intake related constraints for each non-core reaction
        """
        for rxn in self._rncore:
            activation_var = self._activation_vars[rxn]

            # variable that should be bounded by carbon_uptake
            reac_var = rxn.forward_variable + rxn.reverse_variable + activation_var * self._C_uptake
            # adding the constraint to the model
            self._tfa_model.add_constraint(kind=MyConstraintClass,
                                           hook=rxn,
                                           expr=reac_var,
                                           ub=self._C_uptake,
                                           queue=True)

        # push constraints in one bulk (faster)
        self._tfa_model._push_queue()
        # refresh constraint fields
        self._tfa_model.repair()



    def _generate_objective(self):
        """
        Generate and add the maximization objective : set as many activation variables as possible to 1 (deactivated)
        """
        # Sum of binary variables to be maximized
        objective_sum = symbol_sum(list(self._activation_vars.values()))
        # Set the sum as the objective function
        self._tfa_model.objective = self._tfa_model.problem.Objective(objective_sum, direction='max')

    def run_optimisation(self):
        self._tfa_model.prepare()

        # Deactivate tfa computation for non-core reactions
        #TODO : check whether it must be done once or at each iteration
        for ncrxn in self._rncore:
            ncrxn.thermo['computed'] = False

        self._tfa_model.convert()

        tfa_solution = self._tfa_model.optimize()
        return tfa_solution

    def lump_reaction(self, bio_rxn):
        """
        :param bio_rxn: The objective biomass reaction to lump
        :return: the lumped reaction
        """

        # Computing TFA solution
        solution = self.run_optimisation()

        # TODO use generators to improve speed
        # TODO maybe use sympy.add
        lumped_core_reactions = sum([rxn * solution.fluxes.get(rxn.id) for rxn in self._rcore])
        lumped_ncore_reactions = sum([rxn * solution.fluxes.get(rxn.id)*self._activation_vars[rxn].variable.primal for rxn in self._rncore])
        lumped_BBB_reactions = sum([rxn * solution.fluxes.get(rxn.id) for rxn in self._rBBB])

        lumped_reaction = sum([lumped_core_reactions, lumped_ncore_reactions, lumped_BBB_reactions])

        return lumped_reaction



