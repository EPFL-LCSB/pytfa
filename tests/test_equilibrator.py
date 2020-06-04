#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the equilibrator integration."""
import pytest
import sys
from pytfa.io.json import json_dumps_model, json_loads_model



# previous version aren't compatible for the required equlibrator_cache version
# this annotation is not necessary since this file shoud be ignored in the test
# collection phase for previous versions, but it makes it more clear.
@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_load_with_equi_thermo(request):
    """Build thermo data structure with equilibrator."""
    t_model = request.config.cache.get('model', None)
    from pytfa.thermo.equilibrator import build_thermo_from_equilibrator
    from pytfa import ThermoModel
    from settings import cobra_model
    cmodel = cobra_model.copy()
    for met in cmodel.metabolites:
        # normalize but don't overwrite
        met.annotation["seed.compound"] = met.annotation["seed_id"]
    thermo_data = build_thermo_from_equilibrator(cmodel)
    t_model = ThermoModel(thermo_data, cmodel)
    t_model.solver = "optlang-glpk"
    request.config.cache.set('model', json_dumps_model(t_model))



@pytest.mark.dependency(depends=['test_load_with_equi_thermo'])
@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_equilibrator_preparation(request):
    """Prepare using equilibrator_api."""
    t_model = json_loads_model(request.config.cache.get('model',None))
    assert t_model is not None
    t_model.prepare()
    request.config.cache.set('model', json_dumps_model(t_model))


@pytest.mark.dependency(depends=['test_equilibrator_preparation'])
@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_equilibrator_conversion(request):
    """Convert by the usual method in `ThermoModel`."""
    t_model = json_loads_model(request.config.cache.get('model',None))
    assert t_model is not None
    t_model.convert()
    request.config.cache.set('model', json_dumps_model(t_model))


@pytest.mark.dependency(depends=['test_equilibrator_conversion'])
@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_equilibrator_optim(request):
    """LP optimization by the usual method in `ThermoModel`."""
    t_model = json_loads_model(request.config.cache.get('model',None))
    assert t_model is not None
    solution = t_model.optimize()
    return sum(solution.fluxes)
