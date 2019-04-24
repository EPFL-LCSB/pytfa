#!/usr/bin/env python
# -*- coding: utf-8 -*-

from cobra import Reaction
from ..optim.utils import symbol_sum

from pytfa.optim.variables import ReactionVariable, BinaryVariable, get_binary_type
from pytfa.optim.constraints import ReactionConstraint

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
    def __init__(self, tfa_model, biomass_rxns, core_subsystems, carbon_uptake, growth_rate):
        """
        : param tfa_model: The GEM (associated with the thermodynamics constraints) that lumpGEM must work on
        : type GEM: pytfa model

        : param biomass_rxns: list of biomass reactions
        : type biomass_rxns: [GEM.biomass_rxn.id]

        : param core_subsystems: list of Core subsystems names
        : type core_subsystems: [string]

        : param carbon_intake: the amount of carbon atoms the cell intakes from its surrounding
        : type carbon_intake: float

        : param growth_rate: theoretical maximum specific growth rate in 1/hr units
        : type growth_rate: float
        """

        self._tfa_model = tfa_model

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
        # TODO default : solver du modele
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
        self._sinks = self._prepare_sinks()

    def _generate_carbon_constraints(self):
        """
        Generate carbon intake related constraints for each non-core reaction
        For each reaction rxn : rxn.forward_variable + rxn.reverse_variable + activation_var * C_uptake < C_uptake
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

    def _prepare_sinks(self):
        """
        For each BBB (reactant of the biomass reactions), generate a sink, i.e an unbalanced reaction BBB ->
        of which purpose is to enable the BBB to be output of the GEM
        :return: the dict {BBB: sink} containing every BBB (keys) and their associated sinks
        """
        all_sinks = {}
        print("Preparing sinks...")

        for bio_rxn in self._rBBB:
            print(bio_rxn.id)
            for met, stoech_coeff in bio_rxn.metabolites.items():

                # stoech_coeff < 0 indicates that the metabolite is a reactant
                if (stoech_coeff < 0) and (met not in all_sinks.keys()):
                    print("   " + met.id)
                    sink = Reaction("Sink_" + bio_rxn.id + "_" + met.id)
                    sink.name = "Sink_" + bio_rxn.name + "_" + met.name
                    # Subsystem specific to BBB sinks
                    sink.subsystem = "Demand"

                    # A sink is simply a reaction which consumes the BBB
                    sink.add_metabolites({met: -1})
                    # The sinks will be activated later (cf compute_lumps), one at a time
                    sink.knock_out()

                    # The stoechiometric coefficients will be used to define the lower bound of the sink,
                    # thus it must be stored
                    all_sinks[met] = (sink.id, -stoech_coeff)
                    self._tfa_model.add_reactions([sink])

                # reactant already seen
                elif stoech_coeff < 0:
                    # The BBB has already been associated to a sink, so we simply increase the bound of the sink
                    all_sinks[met][1] -= stoech_coeff
                    # Equivalent to this, but there is a knockout :
                    #self._tfa_model.reactions.get_by_id(sinks[met]).lower_bound += self._growth_rate * stoech_coeff

        self._tfa_model.prepare()
        for ncrxn in self._rncore:
            ncrxn.thermo['computed'] = False

        return all_sinks

    def _generate_objective(self):
        """
        Generate and add the maximization objective : set as many activation variables as possible to 1
        When an activation variable is set to 1, the corresponding non-core reaction is deactivated
        """
        # Sum of binary variables to be maximized
        objective_sum = symbol_sum(list(self._activation_vars.values()))
        # Set the sum as the objective function
        self._tfa_model.objective = self._tfa_model.problem.Objective(objective_sum, direction='max')

    def compute_lumps(self):
        """
        For each BBB (reactant of the biomass reaction), add the corresponding sink to the model, then optimize and
        lump the result into one lumped reaction
        :return: The dict {BBB: lump} containing every lumped reactions, associated to their BBBs
        """

        # Must be called before optimization
        self._tfa_model.convert()

        # dict: {metabolite: lumped_reaction}
        lumps = {}

        for met_BBB, (sink_id, stoech_coeff) in self._sinks.items():
            print("Considering :" + met_BBB.id)

            # Activate reaction by setting its lower bound
            self._tfa_model.reactions.get_by_id(sink_id).lower_bound = self._growth_rate * stoech_coeff

            # TODO timeout
            tfa_solution = self._tfa_model.optimize()

            # TODO maybe use sympy.add
            lumped_core_reactions  = sum([rxn * tfa_solution.fluxes.get(rxn.id) for rxn in self._rcore])
            lumped_ncore_reactions = sum([rxn * tfa_solution.fluxes.get(rxn.id) * self._activation_vars[rxn].variable.primal for rxn in self._rncore])
            lumped_BBB_reactions   = sum([rxn * tfa_solution.fluxes.get(rxn.id) for rxn in self._rBBB])

            lumped_reaction = sum([lumped_core_reactions, lumped_ncore_reactions, lumped_BBB_reactions])

            lumps[met_BBB] = lumped_reaction

            # Deactivating reaction by setting both bounds to 0
            self._tfa_model.reactions.get_by_id(sink_id).knock_out()

        return lumps

