# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamic constraints for Flux-Based Analysis of reactions

.. moduleauthor:: pyTFA team

Unit tests for the core functionalities of pytfa

"""

from cobra.test import create_test_model
import os
import pytfa
import pytfa.io

cobra_model = create_test_model("salmonella")
this_directory = os.path.dirname(os.path.realpath(__file__))

# Load the thermo database
thermo_data = pytfa.io.load_thermoDB(this_directory + '/../data/thermo_data.thermodb')
tmodel = pytfa.ThermoModel(thermo_data, cobra_model)



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

    