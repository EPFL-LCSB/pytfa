# -*- coding: utf-8 -*-
"""Thermodynamic information for metabolites from eQuilibrator.

.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team
"""

import equilibrator_cache.compatibility as compat

from equilibrator_cache import create_compound_cache_from_quilt
from .std import TEMPERATURE_0
from ..utils.numerics import BIGM_DG
from ..utils.logger import get_bistream_logger


ccache = None
cc = None
PhasedReaction = None
logger = get_bistream_logger("eQuilibrator_formation")


def build_thermo_from_equilibrator(model, T=TEMPERATURE_0):
    """Build `thermo_data` structure from a cobra Model.

    The structure of the returned dictionary is specified in the pyTFA
    [documentation](https://pytfa.readthedocs.io/en/latest/thermoDB.html).

    :param model: cobra.Model
    :return thermo_data: dict
        to be passed as argument to initialize a `ThermoModel`.

    """
    global ccache
    global cc
    global PhasedReaction
    # lazy loading
    if ccache is None:
        from equilibrator_api import ComponentContribution, Q_
        from equilibrator_api.phased_reaction import PhasedReaction

        ccache = create_compound_cache_from_quilt()
        cc = ComponentContribution(temperature=Q_(str(T) + "K"))
        logger.debug("eQuilibrator loaded.")

    thermo_data = {"name": "eQuilibrator", "units": "kJ/mol", "cues": {}}
    met_to_comps = compat.map_cobra_metabolites(ccache, model.metabolites)
    thermo_data["metabolites"] = [
        compound_to_entry(met, cc) for met in met_to_comps
    ]
    return thermo_data


def compute_dGf(compound, cc):
    """Get ΔGf from equilibrator `compound`."""
    dG0_prime, dG0_uncertainty = cc.dG0_prime(
        PhasedReaction(sparse={compound: 1}, rid=f"tmp_{compound.id}")
    )
    return dG0_prime, dG0_uncertainty


def compound_to_entry(compound, cc):
    """Build thermo structure entry from a `equilibrator_cache.Compound`.

    eQuilibrator works with Component Contribution instead of groups, so it is
    not possible to generate cues from it.

    :param compound: equilibrator_cache.Compound
    :return: dict
        with keys ['deltaGf_std', 'deltaGf_err', 'error', 'struct_cues',
        'id', 'pKa', 'mass_std', 'charge_std', 'nH_std', 'name', 'formula',
        'other_names']

    """
    deltaGf_std, deltaGf_err = (BIGM_DG, BIGM_DG)
    nH_std = compound.atom_bag["H"] if "H" in compound.atom_bag else 0
    try:
        deltaGf_std, deltaGf_err = compute_dGf(compound, cc)
        err = "Nil"
    except Exception as e:
        err = 1
        logger.debug(
            "{} : thermo data NOT created, compound : {}, e : {}".format(
                compound.id, e
            )
        )
    return dict(
        deltaGf_std=deltaGf_std,
        deltaGf_err=deltaGf_err,
        error=err,
        struct_cues=None,
        id=compound.id,
        pKa=compound.dissociation_constants,
        mass_std=compound.mass,
        charge_std=None,
        nH_std=nH_std,
        name=compound.id,
        formula=compound.formula,
        other_names=compound.identifiers,
    )
