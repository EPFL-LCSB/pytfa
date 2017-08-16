#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytfa

from optlang.exceptions import SolverError
from cobra.core.model import SolverNotFound
from cobra.flux_analysis import flux_variability_analysis

from pytfa.io import import_matlab_model, load_thermoDB
from pytfa.optim.relaxation import relax_dgo


CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'


cobra_model = import_matlab_model('../models/small_yeast.mat')


# Load reaction DB
print("Loading thermo data...")

thermo_data = load_thermoDB('../data/thermo_data.thermodb')

print("Done !")

mytfa = pytfa.ThermoModel(thermo_data, cobra_model)
mytfa.normalize_reactions()
mytfa.solver = CPLEX

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


## Try to find an intermediate value
mid_value = (fba_value + tfa_value)/2

mytfa.reactions.biomass_SC5_notrace.lower_bound = mid_value

relaxed_mytfa, slacked_tfa, relax_table = relax_dgo(mytfa, solver = CPLEX)

solver_results = dict()

for solver in [GLPK,CPLEX,GUROBI]:
    try:
        relaxed_mytfa.solver = solver
        this_sol = relaxed_mytfa.optimize()
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
