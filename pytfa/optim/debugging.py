# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Debugging of models

"""

from collections import defaultdict
import pandas as pd

def debug_iis(model):
    """
    Performs reduction to an Irreducible Inconsistent Subsystem (IIS)

    :param model:
    :return:
    """
    out_c = list()
    out_v = list()
    if model.solver.__class__.__module__ == 'optlang.gurobi_interface':
        model.solver.problem.computeIIS()
        print('# Constraints:')
        for c in model.solver.problem.getConstrs():
            if c.IISConstr:
                out_c.append(c)
                print('%s' % c.constrName)
        print('# Variables:')
        for v in model.solver.problem.getVars():
            if v.IISUB != 0 or v.IISLB != 0:
                out_v.append(v)
                print('{}: IISLB = {}, IISUB = {}, (original bounds {}, {})'\
                      .format(v.VarName, v.IISLB, v.IISUB, v.LB, v.UB))
        # model.write("IIS_debug_{}.ilp".format(model.name))
    elif model.solver.__class__.__module__ == 'optlang.cplex_interface':
        model.solver.problem.conflict.refine(
                model.solver.problem.conflict.all_constraints())
    else:
        model.logger.error('Not implemented for solver {}'.format(model.solver))

    return out_c, out_v

def find_extreme_coeffs(model,n=5):
    max_coeff_dict = defaultdict(int)
    min_coeff_dict = defaultdict(lambda:1000)
    max_cons_dict = dict()
    min_cons_dict = dict()

    for the_cons in model.constraints:
        for the_var, the_coeff in the_cons.expression.as_coefficients_dict().items():
            if abs(the_coeff) > max_coeff_dict[the_var.name]:
                max_coeff_dict[the_var.name] = abs(the_coeff)
                max_cons_dict[the_var.name] = the_cons.name
            if 0 < abs(the_coeff) < min_coeff_dict[the_var.name]:
                min_coeff_dict[the_var.name] = abs(the_coeff)
                min_cons_dict[the_var.name] = the_cons.name

    def prep_result(cons_dict, coeff_dict):
        coeff_data = pd.DataFrame.from_dict(coeff_dict, orient = 'index')
        cons_data = pd.DataFrame.from_dict(cons_dict, orient = 'index')

        res = pd.concat([cons_data, coeff_data], axis = 1)
        res.columns = ['constraint','coeff']
        res.index.name = 'variable'
        return res


    ret1 = prep_result(max_cons_dict, max_coeff_dict)\
        .sort_values('coeff',ascending=False).head(n)
    ret2 = prep_result(min_cons_dict, min_coeff_dict)\
        .sort_values('coeff',ascending=True).head(n)

    return pd.concat([ret1, ret2], axis = 0)

def find_maxed_vars(model, ub = 1000, epsilon = 1e-2):
    ret = [(x.name, x.primal, x.ub)
            for x in model.variables
            if abs(x.ub - x.primal) <epsilon and ub>x.ub >0
                                            and not x.type=='binary']

    return pd.DataFrame(ret)