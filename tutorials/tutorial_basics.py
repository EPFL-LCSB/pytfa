#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytfa

from optlang.exceptions import SolverError

from cobra.core.model import SolverNotFound
from cobra.flux_analysis import flux_variability_analysis
from cobra.io import load_matlab_model, load_json_model


from pytfa.io import import_matlab_model, load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data


CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'
solver = GLPK

case = 'reduced' # 'reduced' or full'

# Load reaction DB
print("Loading thermo data...")

thermo_data = load_thermoDB('../data/thermo_data.thermodb')

print("Done !")

if case == 'reduced':
    cobra_model = import_matlab_model('../models/small_ecoli.mat')
    mytfa = pytfa.ThermoModel(thermo_data, cobra_model)
    biomass_rxn = 'Ec_biomass_iJO1366_WT_53p95M'
elif case == 'full':
    # We import pre-compiled data as it is faster for bigger models
    cobra_model = load_json_model('../models/iJO1366.json')

    lexicon = read_lexicon('../models/iJO1366/lexicon.csv')
    compartment_data = read_compartment_data('../models/iJO1366/compartment_data.json')

    # Initialize the cobra_model
    mytfa = pytfa.ThermoModel(thermo_data, cobra_model)

    # Annotate the cobra_model
    annotate_from_lexicon(mytfa, lexicon)
    apply_compartment_data(mytfa, compartment_data)

    biomass_rxn = 'Ec_biomass_iJO1366_WT_53p95M'

mytfa.name = 'tutorial_basics'
mytfa.solver = solver
mytfa.objective = biomass_rxn

# Solver settings

def apply_solver_settings(model, solver = solver):
    model.solver = solver
    # model.solver.configuration.verbosity = 1
    model.solver.configuration.tolerances.feasibility = 1e-9
    if solver == 'optlang_gurobi':
        model.solver.problem.Params.NumericFocus = 3
    model.solver.configuration.presolve = True

apply_solver_settings(mytfa)


## FBA
fba_solution = cobra_model.optimize()
fba_value = fba_solution.objective_value
# fva = flux_variability_analysis(mytfa)

## TFA conversion
mytfa.prepare()
mytfa.convert()#add_displacement = True)

## Info on the cobra_model
mytfa.print_info()

## Optimality
tfa_solution = mytfa.optimize()
tfa_value = tfa_solution.objective_value

# It might happen that the model is infeasible. In this case, we can relax 
# thermodynamics constraints:

if tfa_value < 0.1:
    from pytfa.optim.relaxation import relax_dgo

    mytfa.reactions.get_by_id(biomass_rxn).lower_bound = 0.5*fba_value
    relaxed_model, slack_model, relax_table = relax_dgo(mytfa)

    original_model, mytfa = mytfa, relaxed_model

    print('Relaxation: ')
    print(relax_table)
    
    tfa_solution = mytfa.optimize()
    tfa_value = tfa_solution.objective_value

# Report
print('FBA Solution found : {0:.5g}'.format(fba_value))
print('TFA Solution found : {0:.5g}'.format(tfa_value))


solver_results = dict()

# Try different solvers
for solver in [GLPK,CPLEX,GUROBI]:
    try:
        apply_solver_settings(mytfa,solver)
        this_sol = mytfa.optimize()
        solver_results[solver] = this_sol
        print ("{}: {}".format(solver, this_sol.objective_value))
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
