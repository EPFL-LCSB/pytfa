#!/usr/bin/env python
# -*- coding: utf-8 -*-

from cobra import Reaction
from ..optim.utils import symbol_sum

from pytfa.optim.variables import ReactionVariable, BinaryVariable, get_binary_type
from pytfa.optim.constraints import ReactionConstraint

from numpy import sum, round

from optlang.interface import INFEASIBLE, TIME_LIMIT, OPTIMAL

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'


class InfeasibleExcept(Exception):
    def __init__(self, status, feasibility):
        self.status = status
        self.feasibility = feasibility


class TimeoutExcept(Exception):
    def __init__(self, time_limit):
        self.time_limit = time_limit


class MyVariableClass(ReactionVariable, BinaryVariable):
    prefix = 'VC_'

    def __init__(self, reaction, **kwargs):
        ReactionVariable.__init__(self, reaction,
                                  type=get_binary_type(),
                                  **kwargs)


# Define a new constraint type:
class MyConstraintClass(ReactionConstraint):
    prefix = 'CC_'


def is_exchange(rxn):
    return len(rxn.metabolites) == 1

def trim_epsilon_mets(reaction, epsilon):
    rm_dict = {x:-v for x,v in reaction.metabolites.items() if abs(v)<=epsilon}
    reaction.add_metabolites(rm_dict)

class LumpGEM:
    """
    A class encapsulating the LumpGEM algorithm
    """
    def __init__(self, tfa_model, additional_core_reactions, params):
        """
        :param tfa_model: The GEM (associated with the thermodynamics constraints) that lumpGEM must work on
        :type tfa_model: pytfa model

        :param biomass_rxns: list of biomass reactions
        :type biomass_rxns: [GEM.biomass_rxn.id]

        :param core_subsystems: list of Core subsystems names
        :type core_subsystems: [string]

        :param carbon_uptake: the amount of carbon atoms the cell can uptake from its surrounding
        :type carbon_uptake: float

        :param growth_rate: theoretical maximum specific growth rate in 1/hr units
        :type growth_rate: float

        :param timeout_limit: the maximum amount of time allowed to compute each optimization. Default is 3600s (1 hour)
        :type timeout_limit: float (seconds)
        """

        self._tfa_model = tfa_model

        self._param_dict = params
        self.init_params()

        # Set containing every BBB reaction
        self._rBBB = list()
        # Set containing every exchange reaction
        self._exchanges = list()
        # Set containing every transport reaction
        self._transports = list()
        # Set containing every core reaction
        self._rcore = list()
        # Set containing every non-core reaction
        self._rncore = list()

        # For each reaction
        for rxn in self._tfa_model.reactions:
            # If it's a BBB reaction
            if rxn.id in self.biomass_rxns:
                self._rBBB.append(rxn)
            # If it is an exchange reaction
            elif is_exchange(rxn):
                self._exchanges.append(rxn)
            # If it is a transport reaction
            elif 0:#ch(rxn):
                self._transports.append(rxn)
            # If it's a core reaction
            elif rxn.subsystem in self.core_subsystems:
                self._rcore.append(rxn)
            # If it is part of the intrasubsystem expansion
            elif rxn.id in additional_core_reactions:
                self._rcore.append(rxn)
            # If it's neither BBB nor core, then it's non-core
            else:
                self._rncore.append(rxn)

        # Carbon uptake
        self._C_uptake = self.carbon_uptake
        # Growth rate
        self._growth_rate = self.growth_rate

        # TODO : solver choice
        # TODO default : solver du modele
        #self._solver = 'optlang-cplex'

        self._tfa_model.solver.configuration.timeout = self.timeout_limit
        print("Timeout limit is {}s".format(self.timeout_limit))

        # lumpgem binary variables to deactivate non-core reactions. The reaction is deactivated when the value of
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

    def init_params(self):
        self.core_subsystems = self._param_dict["core_subsystems"]
        self.extracellular_system = self._param_dict["extracellular_system"]
        self.biomass_rxns = self._param_dict["biomass_rxns"]

        self.carbon_uptake = self._param_dict["carbon_uptake"]
        self.growth_rate = self._param_dict["growth_rate"]

        self.small_metabolites = self._param_dict["small_metabolites"]
        self.cofactor_pairs = self._param_dict["cofactor_pairs"]
        # Flatten cofactor_pairs list
        self.cofactors = [cofactor for pair in self.cofactor_pairs for cofactor in pair]
        self.inorganics = self._param_dict["inorganics"]

        self.timeout_limit = self._param_dict["timeout"]

    def _generate_carbon_constraints(self):
        """
        Generate carbon intake related constraints for each non-core reaction
        For each reaction rxn : rxn.forward_variable + rxn.reverse_variable + activation_var * C_uptake < C_uptake
        """

        bigM = self._C_uptake
        for rxn in self._rncore:
            activation_var = self._activation_vars[rxn]

            # variable that should be bounded by carbon_uptake
            reac_var = rxn.forward_variable + rxn.reverse_variable + activation_var * bigM
            # fu = self._tfa_model.forward_use_variable .get_by_id(rxn.id)
            # bu = self._tfa_model.backward_use_variable.get_by_id(rxn.id)
            # reac_var = fu + bu + activation_var
            # adding the constraint to the model
            self._tfa_model.add_constraint(kind=MyConstraintClass,
                                           hook=rxn,
                                           expr=reac_var,
                                           ub=bigM,#1,#bigM,#self._C_uptake,
                                           lb=0,
                                           queue=True)

        # push constraints in one bulk (faster)
        self._tfa_model._push_queue()
        # refresh constraint fields
        self._tfa_model.repair()

    def get_cofactor_adjusted_stoich(self,rxn):
        stoich_dict = {x.id:v for x,v in rxn.metabolites.items()}

        for a,b in self.cofactor_pairs:
            try:
                na = stoich_dict[a] # looks like -54 atp_c
                nb = stoich_dict[b] # looks like +53 adp_c

                n = na+nb # looks like -1

                if n == 0:
                    n = na
                    self._tfa_model.logger.warn(
                        'Cofactor pair {}/{} is equimolar in reaction {}'
                        .format(a,b,rxn.id))
                elif n > 0:
                    self._tfa_model.logger.warn(
                        'Cofactor pair {}/{} looks inverted in reaction {}'
                        .format(a,b,rxn.id))

                stoich_dict[a] =  n # looks like 1
                stoich_dict[b] = -n # looks like -1
            except KeyError:
                pass
        return stoich_dict


    def _prepare_sinks(self):
        """
        For each BBB (reactant of the biomass reactions), generate a sink, i.e an unbalanced reaction BBB ->
        of which purpose is to enable the BBB to be output of the GEM
        :return: the dict {BBB: sink} containing every BBB (keys) and their associated sinks
        """
        all_sinks = {}
        print("Preparing sinks...")

        for bio_rxn in self._rBBB:
            stoich_dict = self.get_cofactor_adjusted_stoich(bio_rxn)
            for met in bio_rxn.metabolites:
                stoech_coeff = stoich_dict[met.id]
                # stoech_coeff < 0 indicates that the metabolite is a reactant
                if (stoech_coeff < 0) and (met not in all_sinks.keys()):
                    sink = Reaction("Sink_" + bio_rxn.id + "_" + met.id)
                    sink.name = "Sink_" + bio_rxn.name + "_" + met.name
                    # Subsystem specific to BBB sinks
                    sink.subsystem = "Demand"

                    # A sink is simply a reaction which consumes the BBB
                    sink.add_metabolites({met: -1})
                    # The sinks will be activated later (cf compute_lumps), one at a time
                    # sink.knock_out()

                    # The stoechiometric coefficients will be used to define the lower bound of the sink,
                    # thus it must be stored
                    all_sinks[met] = (sink.id, -stoech_coeff)
                    self._tfa_model.add_reactions([sink])

                # reactant already seen
                elif stoech_coeff < 0:
                    # The BBB has already been associated to a sink, so we simply increase the bound of the sink
                    all_sinks[met][1] -= stoech_coeff

        # Must be called before changing the reaction.thermo['computed'] values
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

    def compute_lumps(self, force_solve=False):
        """
        For each BBB (reactant of the biomass reaction), add the corresponding sink to the model, then optimize and
        lump the result into one lumped reaction
        :param force_solve: Indicates whether the computations must continue when one lumping yields a status "infeasible"
        :return: The dict {BBB: lump} containing every lumped reactions, associated to their BBBs
        """

        # Must be called before optimization
        self._tfa_model.convert()
        # self._tfa_model.objective_direction = 'min'

        # dict: {metabolite: lumped_reaction}
        lumps = {}
        epsilon_int = self._tfa_model.solver.configuration.tolerances.integrality
        epsilon_flux = self._tfa_model.solver.configuration.tolerances.feasibility

        self._tfa_model.objective_direction = 'max'

        for met_BBB, (sink_id, stoech_coeff) in self._sinks.items():

            print("Considering: " + met_BBB.id)

            sink = self._tfa_model.reactions.get_by_id(sink_id)
            # Activate reaction by setting its lower bound
            prev_lb = sink.lower_bound
            min_prod = self._growth_rate * stoech_coeff
            sink.lower_bound = min_prod

            n_da = self._tfa_model.slim_optimize()


            try:
                # Timeout reached
                if self._tfa_model.solver.status == TIME_LIMIT:
                    raise TimeoutExcept(self._tfa_model.solver.configuration.timeout)
                # Not optimal status -> infeasible
                elif self._tfa_model.solver.status != OPTIMAL:
                    raise InfeasibleExcept( self._tfa_model.solver.status,
                                            self._tfa_model.solver.configuration.tolerances.feasibility)
            except (TimeoutExcept, InfeasibleExcept) as err:
                # If the user want to continue anyway, suits him
                if force_solve:
                    # Raise a warning
                    continue
                else:
                    raise err

            # print('Produced {}'.format(sink.flux),
            #       'with {0:.0f} reactions deactivated'.format(n_da))

            lump_dict = dict()

            for rxn in self._rncore:
                if self._activation_vars[rxn].variable.primal < epsilon_int:
                    lump_dict[rxn] = rxn.flux

            lumped_reaction = sum(rxn * (flux / min_prod) 
                for rxn, flux in lump_dict.items())

            if not lumped_reaction:
                # No need for lump
                self._tfa_model.logger.info('Metabolite {} is produced in enough '
                                        'quantity by core reactions'.format(met_BBB.id))
                continue

            lumped_reaction.id   = sink.id  .replace('Sink_','LUMP_')
            lumped_reaction.name = sink.name.replace('Sink_','LUMP_')
            lumped_reaction.subnetwork = lump_dict

            trim_epsilon_mets(lumped_reaction, epsilon=epsilon_flux)

            lumps[met_BBB] = lumped_reaction

            # Deactivating reaction by setting both bounds to 0
            sink.lower_bound = prev_lb
            # sink.knock_out()

        return lumps
