"""Tests for the equilibrator integration."""

import os
import sys
import pytest
import pytfa.io


this_directory = os.path.dirname(os.path.realpath(__file__))

# Load the cobra_model
cobra_model = pytfa.io.import_matlab_model \
                                    (this_directory \
                                     + '/../models/small_ecoli.mat')
tmodel = None


def test_load_without_thermo():
    """The integration with equilibrator resulted in changes in __init__.
    
    This changes allows to instantiate the model without `thermo_data` as arg.
    """
    global tmodel
    tmodel = pytfa.ThermoModel(None, cobra_model)

@pytest.mark.xfail(sys.version_info < (3, 6),
                   reason="Waiting for integration of cobra in eQuilibrator.")
def test_equilibrator_preparation():
    """Prepare using equilibrator_api."""
    global tmodel
    tmodel.prepare_equilibrator()

@pytest.mark.xfail(sys.version_info < (3, 6),
                   reason="Waiting for integration of cobra in eQuilibrator.")
def test_equilibrator_conversion():
    """Convert by the usual method in `ThermoModel`."""
    global tmodel
    tmodel.convert()

@pytest.mark.xfail(sys.version_info < (3, 6),
                   reason="Waiting for integration of cobra in eQuilibrator.")
def test_equilibrator_optim():
    """LP optimization by the usual method in `ThermoModel`."""
    global tmodel
    tmodel.optimize()
