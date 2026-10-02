#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tutorial for usage of the ModelSEED Biochemistry database with pyTFA.

Instead of the bundled thermodynamic database, pyTFA can build its structure
from a local clone of the ModelSEEDDatabase
(https://github.com/ModelSEED/ModelSEEDDatabase), using the group contribution
energies of its compounds. This suits ModelSEED models, whose metabolites are
ModelSEED compounds (``cpd00002_c0``) in ``c0``/``e0`` compartments, as well as
models annotated with ModelSEED ids, like the ``seed.compound`` annotations of
BiGG models.

Requirements
------------

A clone of the ModelSEEDDatabase, and modelseedpy, whose configuration
(``[biochem] path``) locates the clone. Otherwise, pass its path as ``db_path``
or set the ``MODELSEED_DB_PATH`` environment variable.

.. code:: shell

    pip install .[modelseed]

"""

import pytfa
from cobra.io import load_model

from pytfa.io import apply_compartment_data
from pytfa.optim.relaxation import relax_dgo
from pytfa.thermo.modelseed import (
    annotate_seed_ids,
    build_compartment_data,
    build_thermo_from_modelseed,
)

GLPK = 'optlang-glpk'

# 1. Load the cobra_model. A ModelSEED model, e.g. built with modelseedpy's
# MSBuilder, is used the same way
cobra_model = load_model('textbook')
biomass_rxn = 'Biomass_Ecoli_core'

# 2. Annotate the metabolites with the seed_id pyTFA looks for, from their ids
# or their annotations
missing = annotate_seed_ids(cobra_model)
print('Metabolites without a ModelSEED id:', missing)

# 3. Build the thermodynamic data. Its energies are paired with the protonation
# states of the model's metabolites, so that prepare() leaves the reactions as
# balanced by the model. keep_unpaired=False leaves out the few compounds that
# cannot be paired, instead of letting prepare() re-balance their reactions
thermo_data = build_thermo_from_modelseed(cobra_model)

mytfa = pytfa.ThermoModel(thermo_data, cobra_model)
mytfa.name = 'tutorial_modelseed'
mytfa.solver = GLPK

# 4. Compartment data, with the defaults of modelseedpy's FullThermoPkg: the
# cytosol at pH 7 and -160 mV, the extracellular space at pH 6.5 and 0 mV.
# Each parameter takes a value for all compartments or a dict, for example
# build_compartment_data(mytfa, ph={'c': 7.5}, c_max={'e': 0.1})
apply_compartment_data(mytfa, build_compartment_data(mytfa))

# 5. TFA conversion
mytfa.prepare()
mytfa.convert()

mytfa.print_info()

fba_value = cobra_model.slim_optimize()
tfa_value = mytfa.slim_optimize()

# 6. Thermodynamics may block growth. In this case, relax the standard Gibbs
# energies of the fewest reactions
if tfa_value < 0.1 * fba_value:
    mytfa.reactions.get_by_id(biomass_rxn).lower_bound = 0.5 * fba_value
    relaxed_model, slack_model, relax_table = relax_dgo(mytfa)

    print('Relaxation: ')
    print(relax_table)

    mytfa = relaxed_model
    tfa_value = mytfa.slim_optimize()

print('FBA Solution found : {0:.5g}'.format(fba_value))
print('TFA Solution found : {0:.5g}'.format(tfa_value))
