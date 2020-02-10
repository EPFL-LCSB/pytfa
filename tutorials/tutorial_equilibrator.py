#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tutorial for usage of equilibrator with pyTFA.

Instead of using a local thermodynamic database, pyTFA can use eQuilibrator
(https://gitlab.com/elad.noor/equilibrator-api) to build this structure and use
it to compute the thermodynamics MILP problem.

Requirements
------------

On `setup.cfg`, you can see that equilibrator-api and equilibrator-cache are
extra dependencies of the project. You can check the required versions on this
file or install them directly with pip from the project directory:

.. code:: shell

    pip install .[equilibrator]

.. warning::
    This dependency requires Python version >= 3.6.

"""

import os
import pytfa
import pytfa.io

from pytfa.thermo.equilibrator import build_thermo_from_equilibrator


this_directory = os.path.dirname(os.path.realpath(__file__))

# 1. Load the cobra_model
cobra_model = pytfa.io.import_matlab_model(
    this_directory + "/../models/small_ecoli.mat"
)
# cobra_model.optimizer = "glpk"

# 2. Normalize annotation ids (but don't overwrite them)
for met in cobra_model.metabolites:
    met.annotation["seed.compound"] = met.annotation["seed_id"]

# 3. Build the thermodb strucutre and initialize the model
# this process can take some time
thermo_data = build_thermo_from_equilibrator(cobra_model)
tmodel = pytfa.ThermoModel(thermo_data, cobra_model)
tmodel.solver = "optlang-glpk"

# 4. Now, the process follows the same steps as usual
tmodel.prepare()
tmodel.convert()

# 5. Finally, we can optimize the simulation and print some summary
solution = tmodel.optimize()
print("Some fluxes")
print(solution.fluxes.head())
print("\nOnjective solution:", str(solution.objective_value))
print("\nTotal sum of fluxes:", str(solution.fluxes.sum()))
