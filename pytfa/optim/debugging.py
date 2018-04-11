# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Debugging of models

"""

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

