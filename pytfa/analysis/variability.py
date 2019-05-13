# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Variability analysis

"""

from copy import deepcopy
from functools import partial
from multiprocessing import cpu_count
from multiprocessing.pool import Pool

import pandas as pd
from cobra.core import Reaction
from optlang.exceptions import SolverError
from optlang.interface import INFEASIBLE
from tqdm import tqdm

from ..optim import DeltaG
from ..optim.constraints import ForbiddenProfile
from ..optim.utils import get_direction_use_variables
from ..optim.variables import ForwardUseVariable
from ..utils.logger import get_bistream_logger

CPU_COUNT = cpu_count()
BEST_THREAD_RATIO = int(CPU_COUNT/(4*2))    # Four proc per MILP instance,
                                            # times two threads.

def find_bidirectional_reactions(va, tolerance = 1e-8):
    """
    Returns the ids of reactions that can both carry net flux in the forward or
    backward direction.

    :param va: A variability analysis, pandas Dataframe like so:
                                          maximum       minimum
                6PGLter             -8.330667e-04 -8.330667e-04
                ABUTt2r              0.000000e+00  0.000000e+00
                ACALDt               0.000000e+00  0.000000e+00

    :return:
    """

    return va[va['minimum']*va['maximum'] < -tolerance]


def find_directionality_profiles(tmodel, bidirectional, max_iter = 1e4,
                                 solver = 'optlang-glpk'):
    """
    Takes a ThermoModel and performs enumeration of the directionality profiles

    :param tmodel:
    :param max_iter:
    :return:
    """

    raise(NotImplementedError)

    this_tmodel = deepcopy(tmodel)
    this_tmodel.solver = solver
    profiles = dict()

    iter_count = 0

    bidirectional_reactions = this_tmodel.reactions.get_by_any(bidirectional)

    while this_tmodel.solver.status != INFEASIBLE and iter_count < max_iter:

        try:
            solution = this_tmodel.optimize()
        except SolverError:
            break

        profiles[iter_count] = solution.x_dict
        if iter_count > 0:
            sse = sum((profiles[iter_count-1] - profiles[iter_count])**2)
        else:
            sse =0

        tmodel.logger.info(str(iter_count) + ' - ' + str(sse))

        # active_use_variables = get_active_use_variables(this_tmodel,solution)
        active_use_variables = get_direction_use_variables(this_tmodel,solution)
        bidirectional_use_variables = [x for x in active_use_variables \
                                       if x.reaction in bidirectional_reactions]
        bool_id = ''.join('1' if isinstance(x,ForwardUseVariable) else '0' \
                for x in bidirectional_use_variables)
        Tracer()()
        # Make the expression to forbid this expression profile to happen again
        # FP_1101: FU_rxn1 + FU_rxn2 + BU_rxn3 + FU_rxn4 <= 4-1 = 3
        expr = sum(bidirectional_use_variables)
        this_tmodel.add_constraint(ForbiddenProfile,
                                   hook = this_tmodel,
                                   expr = expr,
                                   id = str(iter_count) + '_' + bool_id,
                                   lb = 0,
                                   ub = len(bidirectional_use_variables)-1)

        iter_count += 1

    return profiles,this_tmodel



def _bool2str(bool_list):
    """
    turns a list of booleans into a string

    :param bool_list: ex: '[False  True False False  True]'
    :return: '01001'
    """
    return ''.join(['1' if x else '0' for x in bool_list])


def _variability_analysis_element(tmodel, var, sense):
    tmodel.objective = var
    tmodel.objective.direction = sense
    sol = tmodel.slim_optimize()
    return sol



def variability_analysis(tmodel, kind='reactions', proc_num = BEST_THREAD_RATIO):
    """
    Performs variability analysis, gicven a variable type

    :param tmodel:
    :param kind:
    :param proc_num:
    :return:
    """

    objective = tmodel.objective

    # If the kind variable is iterable, we perform variability analysis on each,
    # one at a time
    if hasattr(kind, '__iter__') and not isinstance(kind, str):
        va = {}
        for k in kind:
            va[k] = variability_analysis(tmodel, kind=k, proc_num = proc_num)
        df = pd.concat(va.values())
        return df
    elif kind == Reaction or    \
            (isinstance(kind, str) and kind.lower() in ['reaction','reactions']):
        these_vars = {r.id : r for r in tmodel.reactions}
    else:
        these_vars = tmodel.get_variables_of_type(kind)
        these_vars = {x.name : x.variable for x in these_vars}

    tmodel.logger.info('Beginning variability analysis for variable of type {}'    \
                .format(kind))

    results = {'min':{}, 'max':{}}
    for sense in ['min','max']:
        for k,var in tqdm(these_vars.items(), desc=sense+'imizing'):
            tmodel.logger.debug(sense + '-' + k)
            results[sense][k] = _variability_analysis_element(tmodel,var,sense)


    tmodel.objective = objective
    df = pd.DataFrame(results)
    df.rename(columns={'min':'minimum','max':'maximum'}, inplace = True)
    return df


def parallel_variability_analysis(tmodel, kind='reactions', proc_num = BEST_THREAD_RATIO):
    """
    WIP.

    :param tmodel:
    :param kind:
    :param proc_num:
    :return:
    """

    raise(NotImplementedError)

    objective = tmodel.objective

    if kind == Reaction or kind.lower() in ['reaction','reactions']:
        these_vars = tmodel.reactions
    else:
        these_vars = tmodel.get_variables_of_type(kind)

    func = partial(_variability_analysis_element, tmodel)

    pool = Pool(processes=proc_num)
    async_result = pool.map_async(func, these_vars)
    pool.close()
    pool.join()

    # aggregated_result = pd.DataFrame(async_result.get(),
    #                                  columns = ['minimize','maximize'])

    tmodel.objective = objective
    return async_result


def calculate_dissipation(tmodel,solution=None):
    if solution is None:
        solution = tmodel.solution

    reaction_id  = [x.id for x in tmodel.reactions]
    fluxes = solution.x_dict[reaction_id]

    deltag_var = tmodel.get_variables_of_type(DeltaG)
    deltag = pd.Series({x.id:solution.x_dict[x.name]
                                      for x in deltag_var})
    dissipation = fluxes*deltag

    return dissipation