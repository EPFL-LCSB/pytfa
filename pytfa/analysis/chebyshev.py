# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Variability analysis

"""
from tqdm import tqdm
import pandas as pd
import numpy as np
from numpy.linalg import norm
from optlang.interface import Constraint
from pytfa.optim.variables import GenericVariable,ModelVariable
from pytfa.optim.constraints import GenericConstraint

# from ..optim.variables import GenericVariable,ModelVariable

BIGM = 1000

class ChebyshevRadius(ModelVariable):
    """
    Variable representing a Chebyshev radius
    """

    prefix = 'CR_'

def is_inequality(constraint):

    if not isinstance(constraint, Constraint):
        the_cons = constraint.constraint
    else:
        the_cons = constraint

    # If one of the bounds is None, it's an inequality
    return the_cons.lb is None or the_cons.ub is None


def chebyshev_center(model, variables, inplace = False, big_m=BIGM,
                     include = list(), exclude=list()):
    """

    *Warning: Only works with pyTFA variables so far*

    :param model:
    :param variables:
    :param inplace:
    :param big_m:
    :return:
    """

    if not inplace:
        new = model.copy()
        new.optimize()
    else:
        new = model

    vars = get_variables(new, variables)
    include_list = get_cons_var_classes(new, include, type = 'cons')
    exclude_list = get_cons_var_classes(new, exclude, type = 'cons')


    # 0 - Create the Chebyshev radius variable
    r = new.add_variable(kind=ChebyshevRadius,
                         hook=model,
                         id_ ='radius',
                         lb=0,
                         ub=big_m,
                         queue=False)

    # 1 - Find the inequalities associated with the variables
    # of type a_i*x - b_i <= 0

    # Enumerate the constraints, check which ones are:
    #   - Inequalities
    #   - Containing at least 1 of the given variables

    cons_to_edit = dict()

    for cons in tqdm(new._cons_dict.values(), desc='Finding constraints'):
        if type(cons) in exclude_list or type(cons) not in include_list:
            continue
        if not is_inequality(cons):
            continue
        var_intersection = set(cons.expr.free_symbols).intersection(vars)
        if not var_intersection:
            continue

        # 2 - For each inequality, find the norm of the vector of coefficients
        # for the variables
        # ||a_i||_2 = sqrt(sum(x**2 for x in coeffs of variables in this eq))
        a_i = {x:cons.expr.coeff(x) for x in var_intersection}
        a_sq = norm(np.array(list(a_i.values()), dtype=float),ord=2)

        cons_to_edit[cons] = a_sq

    # 3 - Replace the constraint bu the same constraint plus the Chebyshev slack
    # a_i*x - ||a_i||_2 * r - b_i <= 0

    for cons, a_sq in tqdm(cons_to_edit.items(),desc='Editing constraints'):

        new_expr = cons.expr

        if cons.constraint.lb is None:
            # It's a <= 0 constraint
            new_expr += a_sq*r
        elif cons.constraint.ub is None:
            # It's a >=0 constraint
            new_expr -= a_sq*r

        cons.change_expr(new_expr)

    new.logger.info('{} constraints edited'.format(len(cons_to_edit)))
    # 4 - Optimize
    new.repair() # Add the queued constraints
    new.objective.direction = 'max'
    new.objective = r.variable
    new.optimize()

    print('Chebyshev Radius: {}'.format(r.variable.primal))
    if r.variable.ub-r.variable.primal <= new.solver.configuration.tolerances.optimality:
        warn('Chebyshev Radius is close to the upper '
             'bound {}. Change the big_m argument to a bigger one.'.format(r.ub))
    var_values = {k.name:k.primal for k in vars}

    return pd.Series(var_values)


def get_cons_var_classes(model, elements, type):

    if len(elements) == 0:
        return []

    if type.lower().startswith('var'):
        GenericClass = GenericVariable
        model_elements = model._var_kinds
    elif type.lower().startswith('cons'):
        GenericClass = GenericConstraint
        model_elements = model._cons_kinds

    # For safety,
    # Update the variable indices
    model.regenerate_variables()
    model.regenerate_constraints()

    if isinstance(elements[0], str):
        ret = [model_elements[elt][0].__class__ for elt in elements]
    elif isinstance(elements[0], GenericClass):
        ret = elements

    return ret



def get_variables(model, variables):
    if isinstance(variables[0], str):
        # These are var names, we have to retrieve the optlang variables
        vars = [model.variables.get(x) for x in variables]
    elif isinstance(variables[0], GenericVariable):
        # These are pyTFA variables, we have to retrieve the optlang variables
        vars = [x.variable for x in variables]
    return vars

