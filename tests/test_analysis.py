import os

from settings import cobra_model, tmodel

from pytfa.analysis.variability import variability_analysis
from pytfa.analysis.manipulation import apply_reaction_variability, \
    apply_generic_variability, apply_directionality

from cobra.flux_analysis import flux_variability_analysis

solution = tmodel.optimize()

def test_va():
    va = flux_variability_analysis(cobra_model)
    with tmodel as m:
        apply_reaction_variability(m, va)


def test_apply_dir():
    with tmodel as m:
        apply_directionality(m,solution)
        m.optimize()