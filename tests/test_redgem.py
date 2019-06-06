from pytfa.redgem.redgem import RedGEM

from pytfa.io import import_matlab_model
from pytfa.io.base import load_thermoDB
from pytfa.thermo.tmodel import ThermoModel
from pytfa.io import  read_compartment_data, apply_compartment_data, read_lexicon, annotate_from_lexicon
from settings import this_directory
from os.path import join

import pytest

# Check if we are running on Travis CI, to make the run lighter
import os
is_travis = 'TRAVIS' in os.environ \
            or 'CI' in os.environ \
            or os.environ.get('USER')=='travis'
DEBUG = False

if is_travis or DEBUG:
    # # Remove 80% of the mets of the biomass reaction so that less lumps need to be computed:
    # print('Travis env detected. Trimming the biomass reaction')
    # bio_rxn = model.reactions.get_by_id('Ec_biomass_iJO1366_WT_53p95M')
    # bio_rxn.add_metabolites({k:-v for e,(k,v) in
    #                          enumerate(bio_rxn.metabolites.items())
    #                          if e != 1})
    #                          # if e%10 != 0})
    from cobra.test import create_test_model
    model = create_test_model('textbook')
    gly_rxns = ['ENO','FBA','FBP','GAPD','PDH','PFK','PGI','PGK','PGM','PPS',
                'PYK','TPI']
    cofactor_regen_rxns = ['NADTRHD','NADH16','NADTRHD','ATPM','ATPS4r']
    for x in gly_rxns:
        model.reactions.get_by_id(x).subsystem = 'Glycolysis'
    for x in cofactor_regen_rxns:
        model.reactions.get_by_id(x).subsystem = 'Cofactor Regeneration'

    path_to_params = join(this_directory, '..', 'tests/redgem_params_textbook.yaml')


else:
    path_to_model = join(this_directory, '..', 'models/small_ecoli.mat')
    model = import_matlab_model(path_to_model)
    path_to_params = join(this_directory, '..', 'tests/redgem_params.yml')

thermoDB = join(this_directory, '..', 'data/thermo_data.thermodb')
path_to_lexicon = join(this_directory, '..', 'models/iJO1366/lexicon.csv')
path_to_compartment_data = join(this_directory, '..', 'models/iJO1366/compartment_data.json')


# Scaling to avoid numerical errors with bad lumps
for rxn in model.reactions:
    if rxn.id.startswith('LMPD_'):
        rxn.add_metabolites({x:v*(0.1 - 1) for x,v in rxn.metabolites.items()})

thermo_data = load_thermoDB(thermoDB)
lexicon = read_lexicon(path_to_lexicon)
compartment_data = read_compartment_data(path_to_compartment_data)

tfa_model = ThermoModel(thermo_data, model)
annotate_from_lexicon(tfa_model, lexicon)
apply_compartment_data(tfa_model, compartment_data)

tfa_model.name = 'Lumped Model'
tfa_model.prepare()
tfa_model.convert()

# tfa_model.solver.configuration.verbosity = True
tfa_model.logger.setLevel = 30

def test_redgem():
    redgem = RedGEM(tfa_model, path_to_params, False)
    rgem = redgem.run()
    obj_val  = rgem.slim_optimize()
    # assert(obj_val > 0)
    return rgem, redgem


if __name__ == '__main__':
    rgem, redgem = test_redgem()