# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

MILP-fu to reformulate problems

"""
import sympy
# import optlang
from collections import namedtuple

# Faster than optlang Constraint object
ConstraintTuple = namedtuple('ConstraintTuple',['name','expression','ub','lb'])

OPTLANG_BINARY = 'binary'

def subs_bilinear(expr):
    """
    Substitutes bilinear forms from an expression with dedicated variables
    :param expr:
    :return:
    """

    bilinear_ix = [isinstance(x,sympy.Mul) for e,x in enumerate(expr.args)]

    new_expr = expr.copy()

    replacement_dict = dict()

    for bix in bilinear_ix:
        term = expr.args[bix]
        name = '__MUL__'.join(term.args)
        z = sympy.Symbol(name = name)

        new_expr = new_expr.subs(term,z)
        replacement_dict[term] = z

    return new_expr, replacement_dict


def glovers_linearization(b, fy, z = None, L=0, U=1000):
    """
    Glover, Fred.
    "Improved linear integer programming formulations of nonlinear integer problems."
    Management Science 22.4 (1975): 455-460.

    Performs Glovers Linearization of a product
    z = b*f(y) <=> z - b*f(y) = 0
    <=>
    {   L*b <= z <= U*b
    {   f(y) - U*(1-b) <= z <= f(y) - L*(1-b)

    where :
    * b is a binary variable
    * f a linear combination of continuous or integer variables y

    :param b:   Must be a binary optlang variable
    :param z:   Must be an optlang variable. Will be mapped to the product so
                that z = b*f(y)
    :param fy:  Must be an expression or variable
    :param L:   minimal value for fy
    :param U:   maximal value for fy
    :return:
    """

    assert(b.type == OPTLANG_BINARY)

    if z is None:
        name = '__MUL__'.join([b.name, fy.name])
        z = sympy.Symbol(name = name)

    # 1st Glovers constraint
    # L*b <= z
    # 0 <= z - L*b
    cons1 = optlang.Constraint(name = name + '_1',
                               expression = z - L*b,
                               lb = 0)
    # 2nd Glovers constraint
    # z <= U*b
    # 0 <= U*b - z
    cons2 = optlang.Constraint(name = name + '_2',
                               expression = U*b - z,
                               lb=0)

    # 3rd Glovers constraint
    # fy - U*(1-b) <= z
    # 0 <= z - fy + U*(1-b)
    cons3 = optlang.Constraint(name = name + '_3',
                               expression = z - fy + U*(1-b),
                               lb = 0)
    # 4th Glovers constraint
    # z <= fy - L*(1-b)
    # 0 <= fy - L*(1-b) - z
    cons4 = optlang.Constraint(name = name + '_4',
                               expression = fy - L*(1-b) - z,
                               lb=0)

    return z, [cons1,cons2,cons3,cons4]


def petersen_linearization(b, x, z = None, M=1000):
    """
    PETERSEN, C,,
    "A Note on Transforming the Product of Variables to Linear Form in Linear CLIFFORD Programs,"
    Working Paper, Purdue University, 1971.

    Performs Petersen Linearization of a product
    z = b*x <=> z - b*x = 0
    <=>
    {   x + M*b - M <= z <= M*b
    {   z <= x

    where :
    * b is a binary variable
    * f a linear combination of continuous or integer variables y

    :param x:   Must be an expression or variable
    :param b:   Must be a binary optlang variable
    :param z:   Must be an optlang variable. Will be mapped to the product so
                that z = b*f(y)
    :param M:   big-M constraint
    :return:
    """

    assert(b.type == OPTLANG_BINARY)

    if z is None:
        name = '__MUL__'.join([b.name, x.name])
        z = sympy.Symbol(name = name)
    else:
        name = z.name

    # 1st Petersen constraint
    # x + M*b - M <= z
    # x + M*b - z <= M
    cons1 = ConstraintTuple(name = name + '_1',
                            expression = x + M*b - z,
                            lb=0,
                            ub = M)
    # 2nd Petersen constraint
    # z <= M*b
    # 0 <= M*b - z
    cons2 = ConstraintTuple(name = name + '_2',
                            expression = M*b - z,
                            lb=0,
                            ub=None)

    # 3rd Petersen constraint
    # z <= x
    # 0 <= x - z
    cons3 = ConstraintTuple(name = name + '_3',
                            expression = x - z,
                            lb = 0,
                            ub = None,
                            )

    return z, [cons1,cons2,cons3]




