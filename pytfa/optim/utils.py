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
import sympy
from cobra.core.solution import Solution

from numbers import Number

from .constraints import GenericConstraint
from .variables import ForwardUseVariable, BackwardUseVariable
from .variables import GenericVariable

SYMPY_ADD_CHUNKSIZE = 100
INTEGER_VARIABLE_TYPES = ('binary','integer')

def get_all_subclasses(cls):
    """
    Given a variable or constraint class, get all the subclassses
    that inherit from it

    :param cls:
    :return:
    """
    all_subclasses = []

    for subclass in cls.__subclasses__():
        all_subclasses.append(subclass)
        all_subclasses.extend(get_all_subclasses(subclass))

    return all_subclasses

def chunk_sum(variables):
    """
    This functions handles the sum of many sympy variables by chunks, which
    somehow increases the speed of the computation

    You can test it in IPython:
    ```python
    a = sympy.symbols('a0:100')
    %timeit (sum(a))
    # >>> 198 µs ± 11.4 µs per loop (mean ± std. dev. of 7 runs, 1 loop each)

    b = sympy.symbols('b0:1000')
    %timeit (sum(b))
    # >>> 1.85 ms ± 356 µs per loop (mean ± std. dev. of 7 runs, 1 loop each)

    c = sympy.symbols('c0:3000')
    %timeit (sum(c))
    # >>> 5min 7s ± 2.57 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
    ```

    See the `github thread <https://github.com/sympy/sympy/issues/13945>`_

    :param variables:
    :return:
    """

    partial_sums = []

    for chunk_no in range(len(variables) % SYMPY_ADD_CHUNKSIZE):
        first_index = chunk_no * SYMPY_ADD_CHUNKSIZE
        last_index = (chunk_no + 1) * SYMPY_ADD_CHUNKSIZE
        this_chunk = variables[first_index:last_index]
        this_sum = sum(this_chunk)
        partial_sums.append(this_sum)

    return sum(partial_sums)

def symbol_sum(variables):
    """
    
    ``` python
    a = symbols('a0:100')
    
    %timeit Add(*a)
    # >>> 10000 loops, best of 3: 34.1 µs per loop
    
    b = symbols('b0:1000')
    
    %timeit Add(*b)
    # >>> 1000 loops, best of 3: 343 µs per loop
    
    c = symbols('c0:3000')
    
    %timeit Add(*c)
    # >>> 1 loops, best of 3: 1.03 ms per loop
    ```
    
    See the `github thread <https://github.com/sympy/sympy/issues/13945>`_
    :param variables:
    :return:
    """

    from sympy import Add
    
    k=0
    # If we encounter a zero, which is a special type, increase k
    while isinstance(variables[k], sympy.numbers.Zero) and k<len(variables):
        k+=1
        if k == len(variables):
            # everything is 0
            return 0

    if k>len(variables): #it's only zeroes
        return 0

    if isinstance(variables[k], GenericVariable):
        return Add(*[x.variable for x in variables])
    elif isinstance(variables[k], optlang.interface.Variable) or    \
         isinstance(variables[k], sympy.Mul) or \
         isinstance(variables[k], sympy.Add) or \
         isinstance(variables[k], Number):
        return Add(*variables)
    else:
        raise ValueError('Arguments should be of type Number, sympy.Add, or sympy.Mul, '
                         'or optlang.Variable, or GenericVariable')


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
        ret = solution.raw[var_ids]
        ret = ret.index.replace(var2rxn)
        return ret
    else:
        return solution.raw[var_ids]

def compare_solutions(models):
    """
    returns the solution dictionnary for each cobra_model
    :param (iterable (pytfa.thermo.ThermoModel)) models:
    :return:
    """
    return pd.concat([x.solution.raw for x in models], axis=1)

def evaluate_constraint_at_solution(constraint, solution):
    """

    :param expression:
    :param solution: pandas.DataFrame, with index as variable names
    :return:
    """

    if isinstance(solution,Solution):
        solution = solution.raw
    if isinstance(constraint, GenericConstraint):
        constraint = constraint.constraint

    # subs_dict = {x:solution.loc[x.name] for x in constraint.variables}
    # return constraint.expression.subs(subs_dict)

    coefs = constraint.get_linear_coefficients(constraint.expression.free_symbols)
    values = {x:solution.loc[x.name] for x in constraint.expression.free_symbols}

    return symbol_sum([coefs[x]*values[x] for x in coefs])


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

    return [x for x in use_variables if abs(solution.raw[x.name]-1)<epsilon]


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

    return [fwd_use_variables.get_by_id(x.id) if solution.raw[x.id] > epsilon
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
    continuous_model = tmodel.copy()
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

def copy_solver_configuration(source, target):
    """
    Copies the solver configuration from a source model to a target model
    :param source:
    :param target:
    :return:
    """

    # LP method
    try:
        target.solver.configuration.lp_method = source.solver.configuration.lp_method
    except AttributeError:
        pass

    # Presolve
    target.solver.configuration.presolve = source.solver.configuration.presolve

    # Timeout
    target.solver.configuration.timeout = source.solver.configuration.timeout

    # Tolerances
    for tol_name in dir(source.solver.configuration.tolerances):
        tol = getattr(source.solver.configuration.tolerances, tol_name)
        setattr(target.solver.configuration.tolerances, tol_name, tol)

    # Additionnal solver-specific settings
    try:
        # Gurobi
        if source.solver.interface.__name__ == 'optlang.gurobi_interface':
            from gurobipy import GurobiError
            for k in dir(source.solver.problem.Params):
                if not k.startswith('_'):
                    try:
                        v = getattr(source.solver.problem.Params, k)
                        setattr(target.solver.problem.Params, k, v)
                    except GurobiError:
                        pass
    except ModuleNotFoundError:
        pass

    # Verbosity
    target.solver.configuration.verbosity = source.solver.configuration.verbosity


