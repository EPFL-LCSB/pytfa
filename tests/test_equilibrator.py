"""Tests for the equilibrator integration."""
equi = pytest.importorskip(  # noqa: F821
    "pytfa.thermo.equilibrator"
)


import os
import pytest
import pytfa
import pytfa.io
import sys

this_directory = os.path.dirname(os.path.realpath(__file__))

# Load the cobra_model
cobra_model = pytfa.io.import_matlab_model(
    this_directory + "/../models/small_ecoli.mat"
)
cobra_model.optimizer = "glpk"
tmodel = None


# previous version aren't compatible with the equlibrator_cache version needed.
# this annotation is not necessary since the importorskip statement should make
# the tests fail before, but this way it is more explicit
@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_load_with_equi_thermo():
    """Build thermo data structure with equilibrator."""
    global tmodel
    global cobra_model
    for met in cobra_model.metabolites:
        # normalize but don't overwrite
        met.annotation["seed.compound"] = met.annotation["seed_id"]
    thermo_data = equi.build_thermo_from_equilibrator(cobra_model)
    tmodel = pytfa.ThermoModel(thermo_data, cobra_model)
    tmodel.solver = 'optlang-glpk'


@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_equilibrator_preparation():
    """Prepare using equilibrator_api."""
    global tmodel
    tmodel.prepare()


@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_equilibrator_conversion():
    """Convert by the usual method in `ThermoModel`."""
    global tmodel
    tmodel.convert()


@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_equilibrator_optim():
    """LP optimization by the usual method in `ThermoModel`."""
    global tmodel
    solution = tmodel.optimize()
    return sum(solution.fluxes)
