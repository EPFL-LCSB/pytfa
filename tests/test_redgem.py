from pytfa.redgem.redgem import RedGEM

from pytfa.io import import_matlab_model
from pytfa.io.base import load_thermoDB
from pytfa.thermo.tmodel import ThermoModel
from pytfa.io import  read_compartment_data, apply_compartment_data, read_lexicon, annotate_from_lexicon


path_to_model = 'models/small_ecoli.mat'
thermoDB = "data/thermo_data.thermodb"
carbon_uptake = 60000
growth = 0.5

model = import_matlab_model(path_to_model)
thermo_data = load_thermoDB(thermoDB)
lexicon = read_lexicon('models/iJO1366/lexicon.csv')
compartment_data = read_compartment_data('models/iJO1366/compartment_data.json')

tfa_model = ThermoModel(thermo_data, model)
annotate_from_lexicon(tfa_model, lexicon)
apply_compartment_data(tfa_model, compartment_data)

tfa_model.name = 'Lumped Model'

path_to_params = 'tests/redgem_params.yml'

redgem = RedGEM(model, path_to_params, False)
redgem.run()