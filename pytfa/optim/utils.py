# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Relaxation of models with constraint too tight

"""
from copy import deepcopy

import optlang
import pandas as pd
from cobra.core.solution import Solution

from .constraints import GenericConstraint
from .variables import ForwardUseVariable, BackwardUseVariable
from .variables import GenericVariable

INTEGER_VARIABLE_TYPES = ('binary','integer')


def get_solution_value_for_variables(solution, these_vars, index_by_reaction = False):
    if isinstance(these_vars[0],GenericVariable):
        var_ids = [x.name for x in these_vars]
    elif isinstance(these_vars[0],optlang.Variable):
        var_ids = [x.name for x in these_vars]
    elif isinstance(these_vars[0],str):
        var_ids = these_vars
    else:
        raise TypeError('''<these_vars> should be of type pytfa.solving.Variable,
                            optlang.Variable, or str''')


    if index_by_reaction:
        var2rxn = {v.name:v.id for v in these_vars}
        ret = solution.x_dict[var_ids]
        ret = ret.index.replace(var2rxn)
        return ret
    else:
        return solution.x_dict[var_ids]

def compare_solutions(models):
    """
    returns the solution dictionnary for each cobra_model
    :param (iterable (pytfa.core.ThermoModel)) models:
    :return:
    """
    return pd.concat([x.solution.x_dict for x in models], axis=1)

def evaluate_constraint_at_solution(constraint, solution):
    """

    :param expression:
    :param solution: pandas.DataFrame, with index as variable names
    :return:
    """

    if isinstance(solution,Solution):
        solution = solution.x_dict
    if isinstance(constraint, GenericConstraint):
        constraint = constraint.constraint

    subs_dict = {x:solution.loc[x.name] for x in constraint.variables}
    return constraint.expression.subs(subs_dict)

def get_active_use_variables(tmodel,solution):
    """
    Returns the active use variables in a solution. Use Variables are binary
    variables that control the directionality of the reaction
     ex:
     FU_ACALDt
     BU_PFK

    :type tmodel: pytfa.core.ThermoModel
    :param tmodel:
    :param solution:
    :return:
    """
    use_variables = tuple(tmodel.get_variables_of_type(BackwardUseVariable)) \
                  + tuple(tmodel.get_variables_of_type(ForwardUseVariable))

    epsilon = tmodel.solver.configuration.tolerances.integrality

    return [x for x in use_variables if abs(solution.x_dict[x.name]-1)<epsilon]


def get_direction_use_variables(tmodel,solution):
    """
    Returns the active use variables in a solution. Use Variables are binary
    variables that control the directionality of the reaction
    The difference with get_active_use_variables is that variables with both
    UseVariables at 0 will return as going forwards. This is to ensure that the
    output size of the function is equal to the number of FDPs
     ex:
     FU_ACALDt
     BU_PFK

    :type tmodel: pytfa.core.ThermoModel
    :param tmodel:
    :param solution:
    :return:
    """
    fwd_use_variables = tmodel.get_variables_of_type(ForwardUseVariable)
    bwd_use_variables = tmodel.get_variables_of_type(BackwardUseVariable)

    epsilon = tmodel.solver.configuration.tolerances.feasibility

    return [fwd_use_variables.get_by_id(x.id) if solution.x_dict[x.id] > epsilon
            else bwd_use_variables.get_by_id(x.id)
            for x in tmodel.reactions ]

def get_primal(tmodel, vartype, index_by_reactions = False):
    """
    Returns the primal value of the cobra_model for variables of a given type
    :param tmodel:
    :param vartype: Class of variable. Ex: pytfa.optim.variables.ThermoDisplacement
    :param index_by_reactions: Set to true to get reaction names as index instead of
        variables. Useful for Escher export
    :return:
    """

    the_vars = tmodel.get_variables_of_type(vartype)

    if index_by_reactions:
        return pd.Series({x.id:x.variable.primal
                          for x in the_vars})
    else:
        return pd.Series({x.name:x.variable.primal
                          for x in the_vars})


def strip_from_integer_variables(tmodel):
    """
    Removes all integer and binary variables of a cobra_model, to make it sample-able
    :param tmodel:
    :return:
    """
    continuous_model = deepcopy(tmodel)
    continuous_model.name = tmodel.name + ' - continuous'

    integer_variables = set()

    constraints_with_integer_variables = []

    # We go through all the constraint descriptors and check if at least one of
    # their variables is in the integer variable list
    for this_cons in continuous_model._cons_dict.values():
        has_integer_variable = False
        for this_var in this_cons.constraint.variables:
            if this_var.type in INTEGER_VARIABLE_TYPES:
                has_integer_variable += True
                this_var_descriptor = continuous_model._var_dict[this_var.name]
                integer_variables.add(this_var_descriptor)
        constraints_with_integer_variables.append(this_cons)

    for this_cons in constraints_with_integer_variables:
        continuous_model.remove_constraint(this_cons)

    for this_var in integer_variables:
        continuous_model.remove_variable(this_var)

    continuous_model.solver.update()
    # This will update the values =
    print('Is the cobra_model still integer ? {}'     \
          .format(continuous_model.solver.is_integer))

    return continuous_model