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
    else:
        model.logger.error('Not implemented for solver {}'.format(model.solver))

    return out_c, out_v

def find_biggest_coeffs(model,n=1):
    coeff_dict = defaultdict(int)
    cons_dict = dict()

    for the_cons in model.constraints:
        for the_var, the_coeff in the_cons.expression.as_coefficients_dict().items():
            if abs(the_coeff) > coeff_dict[the_var.name]:
                coeff_dict[the_var.name] = abs(the_coeff)
                cons_dict[the_var.name] = the_cons.name

    coeff_data = pd.DataFrame.from_dict(coeff_dict, orient = 'index')
    cons_data = pd.DataFrame.from_dict(cons_dict, orient = 'index')

    res = pd.concat([cons_data, coeff_data], axis = 1)
    res.columns = ['constraint','coeff']
    res.index.name = 'variable'
    res = res[res['coeff'] > 1e-10]

    return res.sort_values('coeff',ascending=False).head(n)
