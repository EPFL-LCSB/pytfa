# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Tests the I/O functionalities of pytfa

"""

import os

from settings import small_model, small_tmodel

#############
#    I/O    #
#############

def test_read_write_mat():
    from pytfa.io.base import write_matlab_model, import_matlab_model
    write_matlab_model(small_tmodel, 'tmp.mat')
    import_matlab_model('tmp.mat')
    os.remove('tmp.mat')

def test_read_write_lexicon():
    from pytfa.io.enrichment import write_lexicon, read_lexicon, \
        annotate_from_lexicon
    fname = 'tmp.lex.csv'
    write_lexicon(small_tmodel,fname)
    lex = read_lexicon(fname)
    annotate_from_lexicon(small_model,lex)
    os.remove(fname)

def test_read_write_compdata():
    from pytfa.io.enrichment import write_compartment_data, read_compartment_data, \
        apply_compartment_data
    fname = 'tmp.comp.json'
    write_compartment_data(small_tmodel,fname)
    cd = read_compartment_data(fname)
    apply_compartment_data(small_model,cd)
    os.remove(fname)

def test_io():
    from pytfa.io.dict import model_to_dict, model_from_dict
    from pytfa.optim.utils import copy_solver_configuration

    dictmodel = model_to_dict(small_tmodel)
    new = model_from_dict(dictmodel)

    copy_solver_configuration(small_tmodel, new)

    sol_orig = small_tmodel.slim_optimize()
    sol_new = new.slim_optimize()

    epsilon = 1e-5 # new.solver.configuration.tolerances.optimality is deprecated in GLPK now, see opencobra/optlang/issues/223

    assert(abs(sol_new - sol_orig) < epsilon)

