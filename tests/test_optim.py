# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Tests of the optimization routines of pytfa

"""

from cobra.test import create_test_model
import os
import sys
import pytfa
import pytfa.io

from settings import tmodel
import pytest


def test_variable_addition():
    global tmodel
    from pytfa.optim.variables import DeltaGstd
    reaction0 = tmodel.reactions[0]
    var0 = tmodel.add_variable(DeltaGstd,reaction0,lb = -1000,ub=1000)
    reaction1 = tmodel.reactions[1]
    var1 = tmodel.add_variable(DeltaGstd,reaction1,lb = -1000,ub=1000)
    tmodel.repair()

    the_name = var1.name

    assert the_name in tmodel.variables
    assert var1.reaction == reaction1
    assert var1.id in getattr(tmodel,var1.__attrname__)

    tmodel.remove_variable(var1)
    tmodel.repair()

    assert the_name not in tmodel.variables
    assert var1.id not in getattr(tmodel, var1.__attrname__)

def test_constraint_addition():
    global tmodel

    from pytfa.optim.constraints import NegativeDeltaG
    reaction0 = tmodel.reactions[0]
    cons0 = tmodel.add_constraint(NegativeDeltaG,reaction0,0,lb = -1000,ub=1000)
    reaction1 = tmodel.reactions[1]
    cons1 = tmodel.add_constraint(NegativeDeltaG,reaction1,0,lb = -1000,ub=1000)
    tmodel.repair()

    the_name = cons1.name

    assert the_name in tmodel.constraints
    assert cons1.reaction == reaction1
    assert cons1.id in getattr(tmodel,cons1.__attrname__)

    tmodel.remove_constraint(cons1)
    tmodel.repair()

    assert the_name not in tmodel.constraints
    assert cons1 not in getattr(tmodel, cons1.__attrname__)

@pytest.mark.xfail(sys.version_info < (3, 6),
                   reason="Container updates behave differently")
def test_relax_dgo():
    global tmodel
    from pytfa.optim.relaxation import relax_dgo

    tmodel.reactions.Ec_biomass_iJO1366_WT_53p95M.lower_bound = 1.5
    tmodel.optimize()
    relax_dgo(tmodel)

@pytest.mark.xfail(sys.version_info < (3, 6),
                   reason="Container updates behave differently")
def test_change_expression():
    global tmodel
    cons = list(tmodel._cons_dict.values())[0]
    cons.change_expr(cons.expr + 2)
    tmodel.optimize()
    