# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Pre-tuned configurations for faster solving

"""

def dg_relax_config(model):
    """

    :param model:
    :return:
    """
    # grbtune output on a hard model :
    #
    # Tested 6992 parameter sets in 46793.78s
    #
    # Baseline parameter set: mean runtime 142.09s
    #
    # Improved parameter set 1 (mean runtime 3.27s):
    #
    # 	NormAdjust 0
    # 	BranchDir 1
    # 	DegenMoves 0
    # 	Heuristics 0
    # 	MIPFocus 1
    # 	Cuts 3
    #
    # Improved parameter set 2 (mean runtime 3.30s):
    #
    # 	NormAdjust 0
    # 	BranchDir 1
    # 	DegenMoves 0
    # 	Heuristics 0.001
    # 	PreSparsify 0
    #
    # Improved parameter set 3 (mean runtime 3.34s):
    #
    # 	NormAdjust 0
    # 	BranchDir 1
    # 	DegenMoves 0
    # 	Heuristics 0.001
    #
    # Improved parameter set 4 (mean runtime 5.22s):
    #
    # 	NormAdjust 1
    # 	BranchDir 1
    # 	DegenMoves 0
    #
    # Improved parameter set 5 (mean runtime 7.18s):
    #
    # 	BranchDir 1
    # 	DegenMoves 0

    if model.solver.interface.__name__ == 'optlang.gurobi_interface':
        model.solver.problem.Params.NormAdjust = 0
        model.solver.problem.Params.BranchDir = 1
        model.solver.problem.Params.DegenMoves = 0
        model.solver.problem.Params.Heuristics = 0.001
        model.solver.problem.Params.Cuts = 3
        model.solver.problem.Params.Presolve = 2
        model.solver.problem.Params.Method = 0

