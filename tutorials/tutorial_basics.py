#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytfa

from optlang.exceptions import SolverError
from cobra.core.model import SolverNotFound
from cobra.flux_analysis import flux_variability_analysis

from pytfa.io import import_matlab_model, load_thermoDB


CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'


cobra_model = import_matlab_model('../models/small_ecoli.mat')


# Load reaction DB
print("Loading thermo data...")

thermo_data = load_thermoDB('../data/thermo_data.thermodb')

print("Done !")

mytfa = pytfa.ThermoModel(thermo_data, cobra_model)
mytfa.solver = GLPK

## FBA
fba_solution = mytfa.optimize()
fba_value = fba_solution.f
fva = flux_variability_analysis(mytfa)

## TFA conversion
mytfa.prepare()
mytfa.convert(add_displacement = True)

## Info on the model
mytfa.print_info()

## Optimality
tfa_solution = mytfa.optimize()
tfa_value = tfa_solution.f

# Report
print('FBA Solution found : {0:.5g}'.format(fba_value))
print('TFA Solution found : {0:.5g}'.format(tfa_value))


solver_results = dict()

# Try different solvers
for solver in [GLPK,CPLEX,GUROBI]:
    try:
        mytfa.solver = solver
        this_sol = mytfa.optimize()
        solver_results[solver] = this_sol
        print ("{}: {}".format(solver, this_sol.f))
    except KeyError:
        print("Solver {} not found".format(solver))
    except SolverError as SE:
        print ("Solver {} returned an error:".format(solver))
        print(SE)
    except SolverNotFound as SNF:
        print ("Solver {} not found:".format(solver))
        print(SNF)
    except Exception as E:
        print("An undocumented error happened:")
        print(E)
