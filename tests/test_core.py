# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Tests for the pytfa module


"""

import pytfa.io
import csv
import lpdiff
import pytest
import os
from pytfa.utils import numerics
from settings import tmodel, this_directory, objective_value

# Minimal relative difference between two values to make a test fail
test_precision = 1 * 10 ** -5

metabolites = []

n_mets = 0
with open(this_directory + '/ref/metData.csv') as csvfile:
    columns = {}
    reference_metabolites = csv.reader(csvfile, delimiter=',')
    for row in reference_metabolites:
        # Title line
        if row[0] == 'id':
            for i,_ in enumerate(row):
                columns[row[i]] = i
            continue

        model_met = tmodel.metabolites.get_by_id(row[columns['id']])
        ref_met = {}

        for item in columns:
            ref_met[item] = row[columns[item]]

        metabolites.append({'cobra_model': model_met, 'ref': ref_met})

        n_mets += 1

def relative_error(a, b):
    # max(1, ...) to avoid an issue with null values
    return abs(a - b)/max(1, min(abs(a), abs(b)))

def test_model_nmets():
    # Make sure we test the same number of metabolites
    global n_mets, tmodel
    assert(n_mets == len(tmodel.metabolites))

@pytest.mark.parametrize("metabolite", metabolites)
def test_metabolites_values(metabolite):
    # global test_test_precision
    # Test each enzyme value
    this_model_met = metabolite['cobra_model']
    this_ref_met = metabolite['ref']

    assert(this_model_met.annotation['seed_id'] == this_ref_met['seedid'])
    assert(this_model_met.compartment == this_ref_met['compartment'])

    for thermoval in ['charge_std', 'mass', 'deltaGf_std', 
                        'deltaGf_err','deltaGf_tr']:
        if relative_error(this_model_met.thermo[thermoval], numerics.BIGM_THERMO) < test_precision:
            pass
        else:
            refval = float(this_ref_met[thermoval].replace(',','.'))
            assert(relative_error(this_model_met.thermo[thermoval], refval) < test_precision)


#############
# REACTIONS #
#############

n_rxns = 0
reactions = []
# Compare reactions
with open(this_directory + '/ref/rxnData.csv') as csvfile:
    columns = {}
    rxns = csv.reader(csvfile, delimiter=',')
    for rxn in rxns:
        # Title line
        if rxn[0] == 'id':
            for i,_ in enumerate(rxn):
                columns[rxn[i]] = i
            continue

        model_rxn = tmodel.reactions.get_by_id(rxn[columns['id']])
        ref_rxn = {}

        for item in columns:
            ref_rxn[item] = rxn[columns[item]]

        reactions.append({'cobra_model': model_rxn, 'ref': ref_rxn})

        n_rxns += 1


def test_model_nrxns():
    # Make sure we test the same number of reactions
    # global n_rxns, tmodel
    assert(n_rxns == len(tmodel.reactions))

@pytest.mark.parametrize("reaction", reactions)
def test_reactions_values(reaction):
    # global test_test_precision
    # Test each enzyme value
    model_rxn = reaction['cobra_model']
    rxn = reaction['ref']

    # Test each enzyme value
    assert(relative_error(model_rxn.lower_bound,
              float(rxn['lower_bound'])) < test_precision)
    assert(relative_error(model_rxn.upper_bound,
              float(rxn['upper_bound'])) < test_precision)
    assert(relative_error(model_rxn.objective_coefficient,
              float(rxn['objective'])) < test_precision)
    assert(model_rxn.thermo['isTrans'] == int(rxn['isTrans']))

    for thermoval in ['computed']:
        refval = float(rxn[thermoval].replace(',','.'))
        assert(relative_error(model_rxn.thermo[thermoval], refval) < test_precision)

    for thermoval in ['deltaGR', 'deltaGRerr']:
        refval = float(rxn[thermoval].replace(',','.'))
        if relative_error(model_rxn.thermo[thermoval], numerics.BIGM_DG) < test_precision:
            pass
        else:
            assert(relative_error(model_rxn.thermo[thermoval], refval) < test_precision)

############
# LP FILES #
############
pytfa.io.writeLP(tmodel, this_directory + '/test.lp')

models = [
          lpdiff.parse_file(this_directory + '/test.lp'),
          lpdiff.parse_file(this_directory + '/ref/reference.lp')
         ]

os.remove(this_directory + '/test.lp')

@pytest.mark.skip(reason="WIP")
def test_lpfiles():
    # global models
    assert(lpdiff.compare(models) < 2)

@pytest.mark.xfail(reason="GLPK can be unreliable")
def test_objective_value():
    # global tmodel, objective_value
    tmodel.optimize()
    assert(relative_error(tmodel.objective.value, objective_value) < test_precision)

