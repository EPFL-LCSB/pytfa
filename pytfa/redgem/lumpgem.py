#!/usr/bin/env python
# -*- coding: utf-8 -*-

from cobra import Reaction

from ..optim.utils import symbol_sum
from ..thermo.utils import is_exchange, check_transport_reaction
from .utils import trim_epsilon_mets, round

from ..optim.variables import ReactionVariable, BinaryVariable, get_binary_type
from ..optim.constraints import ReactionConstraint, ForbiddenProfile

from optlang.interface import INFEASIBLE, TIME_LIMIT, OPTIMAL

from tqdm import tqdm

from collections import defaultdict, namedtuple

import sympy

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'


# Transforms (OnePerBBB --> oneperbbb), (one_per_bbb --> oneperbbb), etc ...
disambiguate = lambda s:s.lower().replace('_','')

Lump = namedtuple('Lump', ['id_', 'metabolites', 'subnetwork', 'gene_reaction_rule'])

class InfeasibleExcept(Exception):
    def __init__(self, status, feasibility):
        self.status = status
        self.feasibility = feasibility


class TimeoutExcept(Exception):
    def __init__(self, time_limit):
        self.time_limit = time_limit


class FluxKO(ReactionVariable, BinaryVariable):
    prefix = 'KO_'

    def __init__(self, reaction, **kwargs):
        ReactionVariable.__init__(self, reaction,
                                  type=get_binary_type(),
                                  **kwargs)

# Define a new constraint type:
class UseOrKOInt(ReactionConstraint):
    prefix = 'UKI_'
# Define a new constraint type:
class UseOrKOFlux(ReactionConstraint):
    prefix = 'UKF_'

def rescale_reactions(rxn_list, factor):
    """
    Rescales the reactions by multiplying stoichiometric coefficients by the factor, and dividing the bounds by the same factore.
    useful with mets that have a low stoichiometry in the biomass reaction, causing low fluxes that might not be detected with the
    integer formulation
    :param rxn_list:
    :param factor:
    :return:
    """

    model = rxn_list[0].model
    model.logger.info('Rescaling {} reactions in the model with a factor of {}'.format(len(rxn_list),factor))

    for this_reaction in rxn_list:
        metabolites = this_reaction.metabolites
        new_metabolites = {k: v * factor \
                           for k, v in metabolites.items()}
        this_reaction.add_metabolites(new_metabolites, combine = False)
        this_reaction.lower_bound /= factor
        this_reaction.upper_bound /= factor

class LumpGEM:
    """
    A class encapsulating the LumpGEM algorithm
    """
    def __init__(self, tfa_model, additional_core_reactions, params, bigM=100):

        """
        :param tfa_model: The GEM (associated with the thermodynamics constraints) that lumpGEM must work on
        :type tfa_model: pytfa model

        :param biomass_rxns: list of biomass reactions
        :type biomass_rxns: [GEM.biomass_rxn.id]

        :param core_subsystems: list of Core subsystems names
        :type core_subsystems: [string]

        :param growth_rate: theoretical maximum specific growth rate in 1/hr units
        :type growth_rate: float

        :param timeout_limit: the maximum amount of time allowed to compute each optimization. Default is 3600s (1 hour)
        :type timeout_limit: float (seconds)
        """

        self._tfa_model = tfa_model
        self.bigM = bigM

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
            # If it's a biomass reaction
            if rxn.id in self.biomass_rxns:
                self._rBBB.append(rxn)
            # Priority to additional core reactions
            elif rxn.id in additional_core_reactions:
                self._rcore.append(rxn)
            # If it is an exchange reaction
            elif is_exchange(rxn):
                self._exchanges.append(rxn)
            # If it is a transport reaction
            elif check_transport_reaction(rxn):
                self._transports.append(rxn)
            # If it's a core reaction
            elif rxn.subsystem in self.core_subsystems:
                self._rcore.append(rxn)
            # If it is part of the intrasubsystem expansion
            # elif rxn.id in core_reactions:
            #     self._rcore.append(rxn)
            # If it's neither BBB nor core, then it's non-core
            else:
                self._rncore.append(rxn)

        self._for_lumps = self._rncore + self._transports + self._exchanges # + self._rcore #Debug last one
        self._lumping_solutions = dict() # Keep track of solutions

        # Growth rate
        self._growth_rate = self.growth_rate

        self._tfa_model.solver.configuration.timeout = self.timeout_limit
        print("Timeout limit is {}s".format(self.timeout_limit))

        try: # Gurobi
            # self._tfa_model.solver.problem.Params.MIPFocus = 3 # for integer emphasis
            # self._tfa_model.solver.problem.Params.Presolve = 2 # for integer emphasis
            if 'mip_gap_abs' in self._param_dict['solver_config'] and self._param_dict['solver_config']['mip_gap_abs']>0:
                self._tfa_model.solver.problem.Params.MIPGapAbs = self._param_dict['solver_config']['mip_gap_abs'] # MIP Gap, absolute
            if 'mip_gap_rel' in self._param_dict['solver_config'] and self._param_dict['solver_config']['mip_gap_rel']>0:
                self._tfa_model.solver.problem.Params.MIPGapRel = self._param_dict['solver_config']['mip_gap_rel'] # MIP Gap, relative

        except AttributeError:
            pass

        # lumpgem binary variables to deactivate non-core reactions. The reaction is deactivated when the value of
        # the variable is 1
        self._activation_vars = {rxn: self._tfa_model.add_variable(kind=FluxKO,
                                                                   hook=rxn,
                                                                   lb=0,
                                                                   ub=1,
                                                                   queue=False)
                                 for rxn in self._for_lumps}
                                 # for rxn in self._rncore}

        self._generate_usage_constraints()
        self._generate_objective()
        self._sinks = self._prepare_sinks()

    def init_params(self):
        self.core_subsystems = self._param_dict["core_subsystems"]
        self.extracellular_system = self._param_dict["extracellular_system"]
        self.biomass_rxns = self._param_dict["biomass_rxns"]

        self.growth_rate = self._param_dict["growth_rate"]

        self.small_metabolites = self._param_dict["small_metabolites"]
        self.cofactor_pairs = self._param_dict["cofactor_pairs"]
        # Flatten cofactor_pairs list
        self.cofactors = [cofactor for pair in self.cofactor_pairs for cofactor in pair]
        self.inorganics = self._param_dict["inorganics"]

        self.timeout_limit = self._param_dict["timeout"]

        self.constraint_method = self._param_dict["constraint_method"]

    def _generate_usage_constraints(self):
        """
        Generate carbon intake related constraints for each non-core reaction
        For each reaction rxn : rxn.forward_variable + rxn.reverse_variable + activation_var * C_uptake < C_uptake
        """
        flux_methods = ['flux', 'fluxes', 'both']
        int_methods = ['int', 'integer', 'both']

        epsilon = self._tfa_model.solver.configuration.tolerances.feasibility

        if self.constraint_method.lower() not in flux_methods + int_methods:
            raise ArgumentError('{} is not a correct constraint method. '
                                'Choose among [Flux, Integer, Both]. '
                                'If you do not know what to choose, go for Flux.'
                                'If it is too slow, go for integer.'
                                'If you get strange lumps, go for both'
                                .format(self.constraint_method))

        # for rxn in self._rncore:
        for rxn in self._for_lumps:
            activation_var = self._activation_vars[rxn]
            if self.constraint_method.lower() in flux_methods:
                reac_var = rxn.forward_variable + rxn.reverse_variable + activation_var * self.bigM
                # adding the constraint to the model
                self._tfa_model.add_constraint(kind=UseOrKOFlux,
                                               hook=rxn,
                                               expr=reac_var,
                                               ub=self.bigM,
                                               lb=0,
                                               queue=True)
            if self.constraint_method.lower() in int_methods:
                fu = self._tfa_model.forward_use_variable.get_by_id(rxn.id)
                bu = self._tfa_model.backward_use_variable.get_by_id(rxn.id)
                reac_var = fu + bu + activation_var
                # adding the constraint to the model
                self._tfa_model.add_constraint(kind=UseOrKOInt,
                                               hook=rxn,
                                               expr=reac_var,
                                               ub=1,
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

                if na == 0 or nb == 0:
                    continue # no need to do anything
                elif na == nb:
                    self._tfa_model.logger.warn(
                        'Cofactor pair {}/{} is equimolar in reaction {}'
                            .format(a,b,rxn.id))
                    continue # no need to do anything
                elif na*nb < 0:
                    # One is a reactant, the other is a product
                    n = na+nb # looks like -1
                else: #na*nb > 0
                    # both are reactant, or product, so they are not cofactoring together
                    continue

                if n > 0:
                    n = -n
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
        For each BBB (reactant of the biomass reactions), generate a sink, i.e an unbalanced reaction BBB -> {}
        of which purpose is to enable the BBB to be output of the GEM
        :return: the dict {BBB: sink} containing every BBB (keys) and their associated sinks
        """
        all_sinks = {}
        print("Preparing sinks...")

        for bio_rxn in self._rBBB:
            stoich_dict = self.get_cofactor_adjusted_stoich(bio_rxn)
            for met in bio_rxn.reactants:
                stoech_coeff = stoich_dict[met.id]
                # stoech_coeff < 0 indicates that the metabolite is a reactant
                if (stoech_coeff < 0) and (met not in all_sinks.keys()):
                    sink = Reaction("Sink_" + bio_rxn.id + "_" + met.id)
                    sink.bounds = (0,self.bigM)
                    sink.name = "Sink_" + bio_rxn.name + "_" + met.name
                    # Subsystem specific to BBB sinks
                    sink.subsystem = "Demand"

                    # A sink is simply a reaction which consumes the BBB
                    # sink.add_metabolites({met: -1})
                    sink.add_metabolites({met: stoech_coeff}) # stoech_coeff is < 0

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
        # for ncrxn in self._rncore:
        #     ncrxn.thermo['computed'] = False

        return all_sinks

    def _generate_objective(self):
        """
        Generate and add the maximization objective : set as many activation variables as possible to 1
        When an activation variable is set to 1, the corresponding non-core reaction is deactivated
        """
        # Sum of binary variables to be maximized
        objective_sum = symbol_sum([v for x,v in self._activation_vars.items() if x in self._for_lumps])
        # Set the sum as the objective function
        self._tfa_model.objective = self._tfa_model.problem.Objective(objective_sum, direction='max')

    def compute_lumps(self, force_solve=False, method='OnePerBBB'):
        """
        For each BBB (reactant of the biomass reaction), add the corresponding sink to the model, then optimize and
        lump the result into one lumped reaction
        :param force_solve: Indicates whether the computations must continue when one lumping yields a status "infeasible"
        :return: The dict {BBB: lump} containing every lumped reactions, associated to their BBBs
        """

        # Must be called before optimization
        self._tfa_model.convert()
        # self._tfa_model.objective_direction = 'min'

        scaling_threshold = 1e-2

        # Set the biomass reactions to 0
        for r in self.biomass_rxns:
            self._tfa_model.reactions.get_by_id(r).lower_bound = 0

        epsilon = self._tfa_model.solver.configuration.tolerances.feasibility

        the_method = disambiguate(method)
        print('Lumping method detected: {}'.format(the_method))

        # dict: {metabolite: lumped_reaction}
        lumps = {}

        self._tfa_model.objective_direction = 'max'

        sink_iter = tqdm(self._sinks.items(), desc = 'met')

        for met_BBB, (sink_id, stoech_coeff) in sink_iter:

            # Cute stuff
            sink_iter.set_description('met={}'.format(met_BBB.id[:10]))
            sink_iter.refresh()

            sink = self._tfa_model.reactions.get_by_id(sink_id)

            # Activate reaction by setting its lower bound
            prev_lb = sink.lower_bound
            # prev_ub = sink.upper_bound

            # d[A]/dt = 0 = ... - n_sink * v_sink = ... - n_bio*v_bio
            # v_sink = n_bio/n_sink * v_bio (/!\ here stoech is >0 but n_sink is <0)
            min_prod = self._growth_rate * -1 * stoech_coeff/sink.metabolites[met_BBB] # this is so that, if the Sink stoich changes, the value stays the same.
            sink.lower_bound = min_prod - epsilon
            # sink.upper_bound = 100
            # sink.lower_bound = 1*self._growth_rate - epsilon

            # Model rescale if min prod is low
            # if min_prod < scaling_threshold:
            #     print('Production requirements of {} are low. scaling the model'.format(met_BBB.id))
            #     rescale_reactions(self._tfa_model.reactions,1e-1/min_prod)

            if the_method == 'oneperbbb':
                this_lump = self._lump_one_per_bbb(met_BBB, sink, force_solve)
                lumped_reactions = [this_lump] if this_lump is not None else list()
            elif the_method.startswith('min+'):
                try:
                    p = int(the_method.replace('min+',''))
                except ValueError:
                    raise ValueError('Min+p method must have p as an integer')
                lumped_reactions = self._lump_min_plus_p(met_BBB, sink, p, force_solve)
            elif the_method.startswith('min'):
                lumped_reactions = self._lump_min_plus_p(met_BBB, sink, 0, force_solve)
            else:
                raise ValueError('Lumping method not recognized: {}. '
                                 'Valid methods are '
                                 'OnePerBBB, Min, Min+p, p natural integer'
                                 .format(the_method))

            sink.lower_bound = prev_lb
            # sink.upper_bound = prev_ub

            # De-scale the model
            # if min_prod < scaling_threshold:
            #     rescale_reactions(self._tfa_model.reactions,1e1*min_prod)

            if not lumped_reactions:
                continue

            lumps[met_BBB] = lumped_reactions

            # Deactivating reaction by setting both bounds to 0
            # sink.knock_out()

        self.lumps = lumps
        return lumps

    def _lump_one_per_bbb(self, met_BBB, sink, force_solve):
        """

        :param met_BBB:
        :param sink:
        :param force_solve:
        :return:
        """

        solution = self._tfa_model.optimize()

        try:
            # Timeout reached
            if self._tfa_model.solver.status == TIME_LIMIT:
                raise TimeoutExcept(self._tfa_model.solver.configuration.timeout)
            # Not optimal status -> infeasible
            elif self._tfa_model.solver.status != OPTIMAL:
                raise InfeasibleExcept(self._tfa_model.solver.status,
                                       self._tfa_model.solver.configuration.tolerances.feasibility)

        except (TimeoutExcept, InfeasibleExcept) as err:
            # If the user want to continue anyway, suits him
            if force_solve:
                # Raise a warning
                return None
            else:
                raise err

        # print('Produced {}'.format(sink.flux),
        #       'with {0:.0f} reactions deactivated'.format(n_da))

        lumped_reaction = self._build_lump(met_BBB, sink, solution)

        return lumped_reaction


    def _lump_min_plus_p(self, met_BBB, sink, p, force_solve):
        """

        :param met_BBB:
        :param sink:
        :param force_solve:
        :return:
        """

        epsilon = self._tfa_model.solver.configuration.tolerances.integrality

        try:
            max_lumps =self._param_dict['max_lumps_per_BBB']
        except KeyError:
            # TODO: Put a warning
            max_lumps=10

        lumps = list()

        with self._tfa_model as model:

            activation_vars = model.get_variables_of_type(FluxKO)
            expr = symbol_sum(activation_vars)

            # Solve a first time, obtain minimal subnet
            model.objective = expr
            model.optimize()
            max_deactivated_rxns = model.objective.value

            # Add constraint forbidding subnets bigger than p

            # The lower bound is the max number of deactivated, minus p
            # Which allows activating the minimal number of reactions, plus p
            lb = max_deactivated_rxns - p
            model.add_constraint(kind=ForbiddenProfile,
                                 hook = model,
                                 id_ = 'MAX_DEACT_{}'.format(met_BBB.id),
                                 expr = expr,
                                 lb = lb,
                                 ub = max_deactivated_rxns,
                                 )

            # n_deactivated_reactions = max_deactivated_rxns

            # Use Zero objective function to enumerate the alternatives
            model.objective = sympy.S.Zero

            # While loop, break on infeasibility
            while len(lumps)<max_lumps:

                try:
                    this_lump = self._lump_one_per_bbb(met_BBB, sink, force_solve)
                except (InfeasibleExcept, TimeoutExcept) as e:
                    if force_solve:
                        pass
                    elif len(lumps) == 0:
                        # No solution AND no lump found
                        raise e

                if model.solver.status != OPTIMAL:
                    break
                elif this_lump is None:
                    # Since the solver is optimal, and we caught optim errors before,
                    # Then the BBB is simply produced in enough quantity by the core
                    break

                lumps.append(this_lump)

                # Add constraint forbidding the previous solution
                is_inactivated = [x for x in activation_vars
                               if abs(x.variable.primal-1) < 2*epsilon]

                # expr = symbol_sum(is_inactivated)
                # model.add_constraint(kind=ForbiddenProfile,
                #                      hook = model,
                #                      id_ = '{}_{}_{}'.format(met_BBB.id,
                #                                              n_deactivated_reactions,
                #                                              len(lumps)),
                #                      expr = expr,
                #                      lb = max_deactivated_rxns-p-1,
                #                      ub = n_deactivated_reactions-1,
                #                      )

                # Add constraint forbidding the previous solution
                is_activated = [x for x in activation_vars
                               if abs(x.variable.primal) < 2*epsilon]

                # print('Active reactions: ', len(is_activated))
                # print('Inactive reactions: ', len(is_inactivated))

                expr = symbol_sum(is_activated)
                model.add_constraint(kind=ForbiddenProfile,
                                     hook = model,
                                     id_ = '{}_{}_{}'.format(met_BBB.id,
                                                             max_deactivated_rxns,
                                                             len(lumps)),
                                     expr = expr,
                                     lb = 1,
                                     ub = max_deactivated_rxns
                                     )
            for c in self._tfa_model.get_condtraints_of_type(ForbiddenProfile):
                self._tfa_model.remove_constraint(c)
            self._tfa_model.repair()

        # TODO: Update of dynamic properties not handled yet
        # upon exiting context manager
        model.repair()
        return lumps


    def _build_lump(self, met_BBB, sink, solution):
        """
        This function uses the current solution of self._tfa_model

        :param met_BBB:
        :param sink:
        :return:
        """

        epsilon_int = self._tfa_model.solver.configuration.tolerances.integrality
        epsilon_flux = self._tfa_model.solver.configuration.tolerances.feasibility

        stoich = sink.metabolites[sink.model.metabolites.get_by_id(met_BBB.id)]
        sigma =  sink.flux * -1 * stoich #1
        lump_dict = dict()

        for rxn in self._rncore:
        # for rxn in self._for_lumps:
            # if self._activation_vars[rxn].variable.primal < 0.5:\
            if abs(solution.fluxes[rxn.id]) > epsilon_flux:\
                    # and abs(rxn.flux) > epsilon_flux:#epsilon_int:
                lump_dict[rxn] = solution.fluxes[rxn.id] / sigma
        # lumped_reaction1 = sum([rxn * (flux / sigma)
        #                       for rxn, flux in lump_dict.items()])

        if not lump_dict:
            # No need for lump
            self._tfa_model.logger.info('Metabolite {} is produced in enough '
                                        'quantity by core reactions'.format(met_BBB.id))
            return None

        lumped_reaction = sum_reactions(lump_dict,
                                        id_=sink.id.replace('Sink_', 'LUMP_'),
                                        epsilon = epsilon_flux,
                                        )
        self._lumping_solutions[lumped_reaction.id_] = solution
        return lumped_reaction


def sum_reactions(rxn_dict, id_ = 'summed_reaction', epsilon = 1e-9):
    """
    Keys are reactions
    Values are their multiplicative coefficient
    """
    stoich = defaultdict(int)

    for rxn,flux in rxn_dict.items():
        if abs(flux) < epsilon:
            continue
        for x, coeff in rxn.metabolites.items():
            # stoich[x.id] += coeff * round(flux, epsilon)
            # stoich[x.id] += round(coeff * flux, epsilon)
            stoich[x.id] += coeff * flux

    gpr = ') and ('.join(x.gene_reaction_rule for x in rxn_dict if x.gene_reaction_rule)

    gpr = ('(' + gpr + ')') if gpr else ''

    # stoich = trim_epsilon_mets(stoich, epsilon=5*epsilon) # ballpark when you sum errors

    # new_rxns = {x.copy():v for x,v in rxn_dict.items()}
    # for x in new_rxns:
    #     x.gene_reaction_rule = ''
    #
    # the_lump = sum([x*v for x,v in new_rxns.items() if not v<epsilon])
    # stoich = trim_epsilon_mets(the_lump.metabolites, epsilon=epsilon)
    # gpr = ''

    new = Lump(id_ = id_,
               metabolites = stoich,
               subnetwork = {x.id:v for x,v in rxn_dict.items()},
               gene_reaction_rule=gpr)

    return new

def print_subnet(model, subnet):
    for r, v in subnet.items():
        print(r, '{:.5g}'.format(v), model.reactions.get_by_id(r).reaction)