from pytfa.redgem.redgem import RedGEM

from pytfa.io import import_matlab_model
from pytfa.io.base import load_thermoDB
from pytfa.thermo.tmodel import ThermoModel
from pytfa.io import  read_compartment_data, apply_compartment_data, read_lexicon, annotate_from_lexicon
from settings import this_directory
from os.path import join


path_to_model = join(this_directory,'..','models/small_ecoli.mat')
# path_to_model = join(this_directory,'..','models/GSmodel_Ecoli.mat')
thermoDB = join(this_directory,'..','data/thermo_data.thermodb')
path_to_lexicon = join(this_directory,'..','models/iJO1366/lexicon.csv')
path_to_compartment_data = join(this_directory,'..','models/iJO1366/compartment_data.json')

model = import_matlab_model(path_to_model)


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

path_to_params = join(this_directory,'..','tests/redgem_params.yml')

redgem = RedGEM(tfa_model, path_to_params, False)
rgem = redgem.run()
