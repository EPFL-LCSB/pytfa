from pytfa.redgem.redgem import RedGEM

from pytfa.io import import_matlab_model
from pytfa.io.base import load_thermoDB
from pytfa.thermo.tmodel import ThermoModel
from pytfa.io import  read_compartment_data, apply_compartment_data, \
    read_lexicon, annotate_from_lexicon
from cobra.io import load_json_model

from os.path import join

# Paths
path_to_model = join('..','models','iJO1366.json')
thermoDB = join('..','data','thermo_data.thermodb')
path_to_lexicon = join('..','models','iJO1366','lexicon.csv')
path_to_compartment_data = join('..','models','iJO1366','compartment_data.json')


# FBA
model = load_json_model(path_to_model)

fba_solution = model.optimize()
fba_value = fba_solution.objective_value

# Thermo prep

thermo_data = load_thermoDB(thermoDB)
lexicon = read_lexicon(path_to_lexicon)
compartment_data = read_compartment_data(path_to_compartment_data)

# Initialize the cobra_model
tfa_model = ThermoModel(thermo_data, model)
tfa_model.name = 'Tutorial'

# Annotate the cobra_model
annotate_from_lexicon(tfa_model, lexicon)
apply_compartment_data(tfa_model, compartment_data)

tfa_model.prepare()
tfa_model.convert()

# tfa_model.solver.configuration.verbosity = True
tfa_model.logger.setLevel = 30

tfa_solution = tfa_model.optimize()
tfa_value = tfa_solution.objective_value

# It might happen that the model is infeasible. In this case, we can relax
# thermodynamics constraints:

if tfa_value < 0.1:
    from pytfa.optim.relaxation import relax_dgo

    biomass_rxn = 'Ec_biomass_iJO1366_WT_53p95M'
    tfa_model.reactions.get_by_id(biomass_rxn).lower_bound = 0.9 * fba_value
    relaxed_model, slack_model, relax_table = relax_dgo(tfa_model, in_place=True)

    original_model, tfa_model = tfa_model, relaxed_model

    print('Relaxation: ')
    print(relax_table)

    tfa_solution = tfa_model.optimize()
    tfa_value = tfa_solution.objective_value

path_to_params = 'tuto_redgem_params.yaml'
redgem = RedGEM(tfa_model, path_to_params, False)
rgem = redgem.run()

redgem_solution = rgem.optimize()
redgem_value = redgem_solution.objective_value

# Report
print('FBA Solution found : {0:.5g}'.format(fba_value))
print('TFA Solution found : {0:.5g}'.format(tfa_value))
print('redGEM Solution found : {0:.5g}'.format(redgem_value))
