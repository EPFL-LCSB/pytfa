#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
.. module:: redgem
   :platform: Unix, Windows
   :synopsis: RedGEM Algorithm

.. moduleauthor:: pyTFA team

Debugging
"""

from cobra import Reaction
from pandas import Series

def make_sink(met, ub=100, lb=0):
    rid = 'sink_' + met.id
    try:
        # if the sink already exists 
        # (happens when debugging several times a model)
        new = met.model.reactions.get_by_id(rid)
    except KeyError:
        new = Reaction(id = rid)
        new.add_metabolites({met:-1})
    
    new.lower_bound = lb
    new.upper_bound = ub

    return new

def add_BBB_sinks(model,biomass_rxn_id, ub=100, lb=0):

    bio_rxn = model.reactions.get_by_id(biomass_rxn_id)

    all_BBBs = bio_rxn.reactants
    all_sinks = list()

    for the_BBB in all_BBBs:
        new = make_sink(the_BBB, ub=ub, lb=lb)
        all_sinks.append(new)

    model.add_reactions([x for x in all_sinks if not x.id in model.reactions])
    return [model.reactions.get_by_id(x.id) for x in all_sinks]

def check_BBB_production(model, biomass_rxn_id, verbose = False):

    all_sinks = add_BBB_sinks(model, biomass_rxn_id, lb = 0)

    prod = dict()

    for the_sink in all_sinks:
        with model:
            model.objective = the_sink.id
            prod[the_sink.id] = model.slim_optimize()

    ret = Series(prod)
    if verbose:
        print(ret)
    return ret

def min_BBB_uptake(model,biomass_rxn_id, min_growth_value, verbose=False):
    

    with model:
        all_sinks = add_BBB_sinks(model, biomass_rxn_id, ub = 0, lb = -100)
        # Uptake is negative
        # Min absolute uptake = Max uptake
        bio_rxn = model.reactions.get_by_id(biomass_rxn_id)
        bio_rxn.lower_bound = min_growth_value
        model.objective = sum( - 1* s.reverse_variable for s in all_sinks)
        model.objective_direction = 'max'
        model.optimize()
        ret = Series({r:r.flux for r in all_sinks})

    if verbose:
        print(ret)
    return ret