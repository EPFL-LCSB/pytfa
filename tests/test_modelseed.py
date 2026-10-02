#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the ModelSEED Biochemistry database integration.

The tests using the database are skipped when no clone of the ModelSEEDDatabase
is found (see :func:`pytfa.thermo.modelseed.find_modelseed_database`).
"""
import pytest
from cobra import Metabolite, Model, Reaction
from cobra.io import load_model

import pytfa
from pytfa.io import apply_compartment_data
from pytfa.io.json import json_dumps_model, json_loads_model
from pytfa.optim.variables import DeltaGstd
from pytfa.thermo.metabolite import MetaboliteThermo
from pytfa.thermo.modelseed import (
    annotate_seed_ids,
    build_compartment_data,
    build_thermo_from_modelseed,
    count_atoms,
    find_modelseed_database,
    get_acidic_pkas,
    get_seed_id_candidates,
    pair_with_species,
)

try:
    find_modelseed_database()
    HAS_DATABASE = True
except FileNotFoundError:
    HAS_DATABASE = False

needs_database = pytest.mark.skipif(
    not HAS_DATABASE, reason="ModelSEEDDatabase clone not found"
)

ATP = {
    "id": "cpd00002",
    "formula": "C10H13N5O13P3",
    "charge_std": -3,
    "nH_std": 13,
    "deltaGf_std": -673.85,
    "deltaGf_err": 6.09,
    "mass_std": 504.0,
    "error": "Nil",
    "pKa": [0.88, 2.73, 3.33, 7.42, 12.6],
    "struct_cues": {},
}


def transformed_energy(entry, ph, ionic_strength):
    return MetaboliteThermo(entry, ph, ionic_strength, thermo_unit="kcal/mol").deltaGf_tr


def test_acidic_pkas():
    marvin = {"pKa": "1:12.60;1:3.33", "pKb": "1:-3.03;1:2.21"}
    assert get_acidic_pkas({"pkas": {"Marvin": marvin}}) == [3.33, 12.6]
    assert get_acidic_pkas({"pkas": None, "pka": "1:14:12.60;1:22:3.33"}) == [3.33, 12.6]
    # NH4+ only has a basic pKa
    assert get_acidic_pkas({"pkas": {"Marvin": {"pKa": "", "pKb": "1:8.86"}}}) == []


def test_count_atoms():
    assert count_atoms("C10H13N5O13P3")["H"] == 13
    assert count_atoms("H")["H"] == 1
    assert count_atoms("HgCl2")["H"] == 0


def test_pair_with_species_keeps_transformed_energy():
    # ATP given for its -2 species instead of its -3 one
    paired = pair_with_species(ATP, "C10H14N5O13P3", -2)
    assert paired["nH_std"] == 14
    assert paired["deltaGf_std"] != ATP["deltaGf_std"]
    for ph, ionic_strength in [(7.0, 0.25), (6.5, 0.0), (8.0, 0.1)]:
        assert transformed_energy(paired, ph, ionic_strength) == pytest.approx(
            transformed_energy(ATP, ph, ionic_strength), abs=1e-9
        )


def test_pair_with_species_refuses_other_compounds():
    assert pair_with_species(ATP, ATP["formula"], ATP["charge_std"]) is ATP
    # One oxygen less
    assert pair_with_species(ATP, "C10H13N5O12P3", -3) is None
    # Protons and charge not moving together
    assert pair_with_species(ATP, "C10H14N5O13P3", -3) is None
    assert pair_with_species(ATP, None, -3) is None


def test_seed_id_candidates():
    assert get_seed_id_candidates(Metabolite("cpd00002_c0")) == ["cpd00002"]
    met = Metabolite("glc__D_e")
    met.annotation = {"seed.compound": ["cpd26821", "cpd00027"], "sbo": "SBO:0000247"}
    assert get_seed_id_candidates(met) == ["cpd26821", "cpd00027"]


def test_compartment_data():
    model = Model()
    model.add_metabolites(
        [Metabolite("cpd00009_c0", compartment="c0"), Metabolite("cpd00009_e0", compartment="e0")]
    )

    data = build_compartment_data(model)
    assert (data["c0"]["pH"], data["e0"]["pH"]) == (7.0, 6.5)
    assert data["c0"]["ionicStr"] == 0.25
    # Transport from the extracellular space into the cytosol
    assert data["e0"]["membranePot"] == {"c0": -160.0, "e0": 0.0}
    assert data["c0"]["membranePot"] == {"c0": 0.0, "e0": 160.0}

    data = build_compartment_data(model, ph={"c": 7.5}, potential=0, c_max={"e0": 0.1})
    assert (data["c0"]["pH"], data["e0"]["pH"]) == (7.5, 6.5)
    assert data["e0"]["membranePot"]["c0"] == 0
    assert (data["c0"]["c_max"], data["e0"]["c_max"]) == (0.02, 0.1)


def test_refuses_convention_b_sources():
    with pytest.raises(NotImplementedError):
        build_thermo_from_modelseed(source="eQuilibrator")


@needs_database
def test_thermo_data():
    thermo_data = build_thermo_from_modelseed()
    assert thermo_data["units"] == "kcal/mol"
    metabolites = thermo_data["metabolites"]

    assert metabolites["cpd00067"]["deltaGf_std"] == -9.5
    assert metabolites["cpd00001"]["formula"] == "H2O"
    # ATP is given for its canonical species, the one of its group contribution
    assert metabolites["cpd00002"]["deltaGf_std"] == -673.85
    # D-glucose-6-phosphate is moved from its -1 group contribution species to
    # its canonical -2 one
    g6p = metabolites["cpd00079"]
    assert (g6p["formula"], g6p["charge_std"], g6p["nH_std"]) == ("C6H11O9P", -2, 11)
    assert g6p["deltaGf_std"] != -430.78


@needs_database
def test_leave_out_unpaired():
    model = load_model("textbook")
    annotate_seed_ids(model)
    species = {}
    for met in model.metabolites:
        if "seed_id" in met.annotation:
            species.setdefault(met.annotation["seed_id"], set()).add(
                (met.formula, met.charge)
            )

    kept = build_thermo_from_modelseed(model)["metabolites"]
    paired = build_thermo_from_modelseed(model, keep_unpaired=False)["metabolites"]
    assert set(paired) < set(kept)
    for seed_id, entry in paired.items():
        assert (entry["formula"], entry["charge_std"]) in species[seed_id]


@needs_database
def test_textbook_model():
    model = load_model("textbook")
    annotate_seed_ids(model)
    assert model.metabolites.g6p_c.annotation["seed_id"] == "cpd00079"
    assert "bigg.metabolite" in model.metabolites.g6p_c.annotation

    tmodel = pytfa.ThermoModel(build_thermo_from_modelseed(model), model)
    tmodel.name = "textbook_modelseed"
    apply_compartment_data(tmodel, build_compartment_data(tmodel))
    tmodel.solver = "optlang-glpk"
    tmodel.prepare()
    tmodel.convert()

    pgi = tmodel.reactions.PGI.thermo
    assert pgi["computed"]
    assert pgi["deltaGR"] == pytest.approx(-0.89, abs=1e-6)
    # Compound errors in quadrature, like the database's rxn00558
    assert pgi["deltaGRerr"] == pytest.approx(3.99, abs=0.01)

    # Every computed reaction has a finite standard deltaG range
    for var in tmodel.get_variables_of_type(DeltaGstd):
        assert var.variable.ub - var.variable.lb < 1e3

    growth = tmodel.slim_optimize()
    assert growth == pytest.approx(model.slim_optimize(), rel=1e-6)
    assert json_loads_model(json_dumps_model(tmodel)).slim_optimize() == pytest.approx(growth)


@needs_database
def test_modelseed_transport():
    """A ModelSEED-style model with c0/e0 compartments and a proton symport"""
    model = Model("modelseed_transport")
    pi_e, pi_c = Metabolite("cpd00009_e0", "HO4P", charge=-2, compartment="e0"), \
                 Metabolite("cpd00009_c0", "HO4P", charge=-2, compartment="c0")
    h_e, h_c = Metabolite("cpd00067_e0", "H", charge=1, compartment="e0"), \
               Metabolite("cpd00067_c0", "H", charge=1, compartment="c0")
    symport = Reaction("rxn05312_c0", lower_bound=-1000, upper_bound=1000)
    symport.add_metabolites({pi_e: -1, h_e: -1, pi_c: 1, h_c: 1})
    exchange = Reaction("EX_cpd00009_e0", lower_bound=-10, upper_bound=1000)
    exchange.add_metabolites({pi_e: -1})
    sink = Reaction("SK_cpd00009_c0", upper_bound=1000)
    sink.add_metabolites({pi_c: -1, h_c: -1})
    proton = Reaction("EX_cpd00067_e0", lower_bound=-1000, upper_bound=1000)
    proton.add_metabolites({h_e: -1})
    model.add_reactions([symport, exchange, sink, proton])
    model.objective = sink

    assert annotate_seed_ids(model) == []
    thermo_data = build_thermo_from_modelseed(model)

    def prepare(**compartment_parameters):
        tmodel = pytfa.ThermoModel(thermo_data, model)
        apply_compartment_data(
            tmodel, build_compartment_data(tmodel, **compartment_parameters)
        )
        tmodel.solver = "optlang-glpk"
        tmodel.prepare()
        return tmodel

    tmodel = prepare()
    thermo = tmodel.reactions.rxn05312_c0.thermo
    assert thermo["computed"] and thermo["isTrans"]
    assert tmodel._transport_compartment == "c0"
    # No proton added: the model is balanced for the species it is given for
    assert tmodel.reactions.rxn05312_c0.metabolites[tmodel.metabolites.cpd00067_c0] == 1

    # Pi(2-) and H+ bring a net -1 charge into the cytosol, at -160 mV
    faraday = 23.061  # kcal/(mol V)
    unpolarized = prepare(potential=0).reactions.rxn05312_c0.thermo
    assert thermo["deltaGR"] - unpolarized["deltaGR"] == pytest.approx(
        faraday * -0.160 * -1
    )

    tmodel.convert()
    assert tmodel.slim_optimize() == pytest.approx(10)
