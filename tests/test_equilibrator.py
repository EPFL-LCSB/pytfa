#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the equilibrator integration."""
import os
import pytest
import pytfa
import pytfa.io
import sys

from pytfa.thermo.equilibrator import build_thermo_from_equilibrator

this_directory = os.path.dirname(os.path.realpath(__file__))

# Load the cobra_model
cobra_model = pytfa.io.import_matlab_model(
    this_directory + "/../models/small_ecoli.mat"
)
cobra_model.optimizer = "glpk"
tmodel = None


# previous version aren't compatible for the required equlibrator_cache version
# this annotation is not necessary since this file shoud be ignored in the test
# collection phase for previous versions, but it makes it more clear.
@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
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
