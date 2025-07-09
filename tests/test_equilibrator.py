#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the equilibrator integration."""
import pytest
import sys
from pytfa.io.json import json_dumps_model, json_loads_model
from pytfa.thermo.equilibrator import build_thermo_from_equilibrator
from pytfa import ThermoModel
from settings import small_model



# previous version aren't compatible for the required equlibrator_cache version
# this annotation is not necessary since this file shoud be ignored in the test
# collection phase for previous versions, but it makes it more clear.
@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_load_with_equi_thermo(request):
    """Build thermo data structure with equilibrator."""
    t_model = request.config.cache.get('model', None)
    cmodel = small_model.copy()
    for met in cmodel.metabolites:
        # normalize but don't overwrite
        if "seed_id" in met.annotation:
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
def test_equilibrator_deltag_loaded(request):
    """Test that deltaG values are actually loaded from equilibrator."""
    # Get the original cobra model and rebuild thermo_data to check it directly
    
    cmodel = small_model.copy()
    for met in cmodel.metabolites:
        # normalize but don't overwrite
        if "seed_id" in met.annotation:
            met.annotation["seed.compound"] = met.annotation["seed_id"]
    
    # Build thermo data from equilibrator
    thermo_data = build_thermo_from_equilibrator(cmodel)
    
    print(f"Thermo data structure: {thermo_data['name']}")
    print(f"Units: {thermo_data['units']}")
    
    # Check metabolites in thermo_data
    metabolites_with_deltag = 0
    sample_mets = []
    
    for met_data in thermo_data['metabolites']:
        if 'DeltaGf_tr' in met_data:
            metabolites_with_deltag += 1
            if len(sample_mets) < 5:  # Collect some samples
                sample_mets.append((met_data['id'], met_data['DeltaGf_tr']))
    
    print(f"Found {metabolites_with_deltag} metabolites with deltaG formation values from equilibrator")
    print("Sample metabolites with deltaG values:")
    for met_id, deltag in sample_mets:
        print(f"  {met_id}: {deltag} {thermo_data['units']}")
    
    # There should be at least some metabolites with deltaG values
    assert metabolites_with_deltag > 0, f"No deltaG formation values found in thermo_data! Expected > 0, got {metabolites_with_deltag}"
    
    return metabolites_with_deltag


@pytest.mark.dependency(depends=['test_equilibrator_deltag_loaded'])
@pytest.mark.skipif(sys.version_info < (3, 6), reason="requires >= python3.6")
def test_equilibrator_optim(request):
    """LP optimization by the usual method in `ThermoModel`."""
    t_model = json_loads_model(request.config.cache.get('model',None))
    assert t_model is not None
    solution = t_model.optimize()
    return sum(solution.fluxes)
