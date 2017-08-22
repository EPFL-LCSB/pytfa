# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamic constraints for Flux-Based Analysis of reactions

.. moduleauthor:: pyTFA team

Tests for the pytfa module


"""
import pytfa
import pytfa.io
import csv
import os
import lpdiff
import pytest

############
# SETTINGS #
############
# Minimal relative difference between two values to make a test fail
Precision = 1 * 10 ** -5
# Objective value of the Matlab solution
Objective_value = 0.263412914022098

######## End of Settings ########
this_directory = os.path.dirname(os.path.realpath(__file__))

# Load the thermo database
thermo_data = pytfa.io.load_thermoDB(this_directory + '/DB_AlbertyUpdate.thermodb')

# Load the model
model = pytfa.io.import_matlab_model(this_directory + '/benchmark_model.mat')

# Make your computations on it
mytfa = pytfa.ThermoModel(thermo_data, model)
mytfa.solver = 'cplex'
mytfa.prepare()
mytfa.convert()

metabolites = []

n_mets = 0
with open(this_directory + '/metData.csv') as csvfile:
    columns = {}
    reference_metabolites = csv.reader(csvfile, delimiter=';')
    for row in reference_metabolites:
        # Title line
        if row[0] == 'id':
            for i,_ in enumerate(row):
                columns[row[i]] = i
            continue

        model_met = mytfa.metabolites.get_by_id(row[columns['id']])
        ref_met = {}

        for item in columns:
            ref_met[item] = row[columns[item]]

        metabolites.append({'model': model_met, 'ref': ref_met})

        n_mets += 1

def relative_error(a, b):
    # max(1, ...) to avoid an issue with null values
    return abs(a - b)/max(1, min(abs(a), abs(b)))

def test_model_nmets():
    # Make sure we test the same number of metabolites
    global n_mets, mytfa
    assert(n_mets == len(mytfa.metabolites))

@pytest.mark.parametrize("metabolite", metabolites)
def test_metabolites_values(metabolite):
    global Precision
    # Test each metabolite value
    this_model_met = metabolite['model']
    this_ref_met = metabolite['ref']

    assert(this_model_met.annotation['SeedID'] == this_ref_met['seedid'])
    assert(this_model_met.compartment == this_ref_met['compartment'])

    for thermoval in ['deltaGf_std', 'deltaGf_err',
                       'charge_std', 'mass', 'deltaGf_tr']:
        refval = float(this_ref_met[thermoval].replace(',','.'))
        assert(relative_error(this_model_met.thermo[thermoval], refval) < Precision)


#############
# REACTIONS #
#############

n_rxns = 0
reactions = []
# Compare reactions
with open(this_directory + '/rxnData.csv') as csvfile:
    columns = {}
    rxns = csv.reader(csvfile, delimiter=';')
    for rxn in rxns:
        # Title line
        if rxn[0] == 'id':
            for i,_ in enumerate(rxn):
                columns[rxn[i]] = i
            continue

        model_rxn = mytfa.reactions.get_by_id(rxn[columns['id']])
        ref_rxn = {}

        for item in columns:
            ref_rxn[item] = rxn[columns[item]]

        reactions.append({'model': model_rxn, 'ref': ref_rxn})

        n_rxns += 1


def test_model_nrxns():
    # Make sure we test the same number of reactions
    global n_rxns, mytfa
    assert(n_rxns == len(mytfa.reactions))

@pytest.mark.parametrize("reaction", reactions)
def test_reactions_values(reaction):
    global Precision
    # Test each metabolite value
    model_met = reaction['model']
    met = reaction['ref']

    assert(relative_error(model_rxn.lower_bound,
              float(rxn[columns['lower_bound']])) < Precision)
    assert(relative_error(model_rxn.upper_bound,
              float(rxn[columns['upper_bound']])) < Precision)
    assert(relative_error(model_rxn.objective_coefficient,
              float(rxn[columns['objective']])) < Precision)
    assert(model_rxn.thermo['isTrans'] == int(rxn[columns['isTrans']]))

    for thermoval in ['computed', 'deltaGR', 'deltaGRerr']:
        refval = float(rxn[columns[thermoval]].replace(',','.'))
        assert(relative_error(model_rxn.thermo[thermoval], refval) < Precision)

############
# LP FILES #
############
pytfa.io.writeLP(mytfa, this_directory + '/test.lp')

models = [
          lpdiff.parse_file(this_directory + '/test.lp'),
          lpdiff.parse_file(this_directory + '/reference.lp')
         ]

#os.remove(this_directory + '/test.lp')

def test_lpfiles():
    global models
    assert(lpdiff.compare(models) < 2)

def test_objective_value():
    global mytfa, Objective_value
    mytfa.optimize()
    assert(relative_error(mytfa.objective.value, Objective_value) < Precision)
