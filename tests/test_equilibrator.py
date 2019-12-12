"""Tests for the equilibrator integration."""

import os
import pytfa
import pytfa.io

from pytfa.thermo.equilibrator import build_thermo_from_equilibrator

this_directory = os.path.dirname(os.path.realpath(__file__))

# Load the cobra_model
cobra_model = pytfa.io.import_matlab_model(
    this_directory + "/../models/small_ecoli.mat"
)
cobra_model.optimizer = "glpk"
tmodel = None


def test_load_with_equi_thermo():
    """Build thermo data structure with equilibrator."""
    global tmodel
    global cobra_model
    for met in cobra_model.metabolites:
        # normalize but don't overwrite
        met.annotation["seed.compound"] = met.annotation["seed_id"]
    thermo_data = build_thermo_from_equilibrator(cobra_model)
    tmodel = pytfa.ThermoModel(thermo_data, cobra_model)
    tmodel.solver = 'optlang-glpk'


def test_equilibrator_preparation():
    """Prepare using equilibrator_api."""
    global tmodel
    tmodel.prepare()


def test_equilibrator_conversion():
    """Convert by the usual method in `ThermoModel`."""
    global tmodel
    tmodel.convert()


def test_equilibrator_optim():
    """LP optimization by the usual method in `ThermoModel`."""
    global tmodel
    solution = tmodel.optimize()
    return sum(solution.fluxes)
