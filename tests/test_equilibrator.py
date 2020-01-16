#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the equilibrator integration."""
import pytest
import sys


t_model = None


# previous version aren't compatible for the required equlibrator_cache version
# this annotation is not necessary since this file shoud be ignored in the test
# collection phase for previous versions, but it makes it more clear.
@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_load_with_equi_thermo():
    """Build thermo data structure with equilibrator."""
    from pytfa.thermo.equilibrator import build_thermo_from_equilibrator
    from pytfa import ThermoModel
    from settings import cobra_model

    global t_model
    cmodel = cobra_model.copy()
    for met in cmodel.metabolites:
        # normalize but don't overwrite
        met.annotation["seed.compound"] = met.annotation["seed_id"]
    thermo_data = build_thermo_from_equilibrator(cmodel)
    t_model = ThermoModel(thermo_data, cmodel)
    t_model.solver = "optlang-glpk"


@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_equilibrator_preparation():
    """Prepare using equilibrator_api."""
    global t_model
    t_model.prepare()


@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_equilibrator_conversion():
    """Convert by the usual method in `ThermoModel`."""
    global t_model
    t_model.convert()


@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_equilibrator_optim():
    """LP optimization by the usual method in `ThermoModel`."""
    global t_model
    solution = t_model.optimize()
    return sum(solution.fluxes)
