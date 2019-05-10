from pytfa.redgem.redgem import RedGEM

from pytfa.io import import_matlab_model
from pytfa.io.base import load_thermoDB
from pytfa.thermo.tmodel import ThermoModel
from pytfa.io import  read_compartment_data, apply_compartment_data, read_lexicon, annotate_from_lexicon
from os.path import join

# Check if we are running on Travis CI, to make the run lighter
import os
is_travis = 'TRAVIS' in os.environ

# Paths
path_to_model = join('..','models','GSmodel_Ecoli.mat')
thermoDB = join('..','data','thermo_data.thermodb')
path_to_lexicon = join('..','models','iJO1366','lexicon.csv')
path_to_compartment_data = join('..','models','iJO1366','compartment_data.json')


# Model prep
model = import_matlab_model(path_to_model)

thermo_data = load_thermoDB(thermoDB)
lexicon = read_lexicon(path_to_lexicon)
compartment_data = read_compartment_data(path_to_compartment_data)

tfa_model = ThermoModel(thermo_data, model)
annotate_from_lexicon(tfa_model, lexicon)
apply_compartment_data(tfa_model, compartment_data)

tfa_model.prepare()
tfa_model.convert()

# tfa_model.solver.configuration.verbosity = True
tfa_model.logger.setLevel = 30

path_to_params = 'tuto_redgem_params.yaml'

if __name__ == '__main__':
    redgem = RedGEM(tfa_model, path_to_params, False)
    rgem = redgem.run()