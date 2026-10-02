# -*- coding: utf-8 -*-
"""Thermodynamic information for metabolites from the ModelSEED Biochemistry
database.

.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Builds the ``thermo_data`` structure of :class:`pytfa.ThermoModel` from a local
clone of the `ModelSEEDDatabase <https://github.com/ModelSEED/ModelSEEDDatabase>`_,
and the compartment data of ModelSEED models (``c0``, ``e0``, ...).

Data policy:

* Energies and errors are ``thermodynamics["Group contribution"]`` of the
  ``Biochemistry/compound_*.json`` records. They are in kcal/mol and in the
  database's Convention A (hydrogens accounted in the compound, H+ at
  -9.5 kcal/mol), the formalism the per-compartment transform of pyTFA expects.
  The other sources ship already transformed to pH 7 (Convention B) and would
  be transformed twice, so they are refused.
* An energy belongs to the species the group contribution was computed for
  (``ModelSEED_GroupContribution.tsv``), which is not always the protonation
  state of the compound in the model. pyTFA transforms every species to the
  least protonated one using the pKas, so the energy is shifted along the same
  pKas to the model's species (else the database's canonical one). The
  transformed energy of formation is unchanged, and :meth:`ThermoModel.prepare`
  does not re-balance the model's reactions with protons. Species pyTFA cannot
  relate keep the group contribution one.
* pKas are Marvin's acidic pKas. Its basic pKas are left out, since pyTFA counts
  every pKa below ``max_ph`` as a deprotonation of the neutral species.
* No structural cues are provided, so :meth:`ThermoModel.prepare` combines the
  compound errors in quadrature, as the database does for its reactions.

modelseedpy provides the default location of the database. Its biochemistry
loader is not used for the values since it drops the per-source
``thermodynamics`` and ``pkas`` fields.
"""

import glob
import json
import logging
import os
import re
from collections import Counter
from copy import deepcopy
from functools import lru_cache

from . import std
from .metabolite import MetaboliteThermo
from ..utils.numerics import BIGM_DG

GROUP_CONTRIBUTION = "Group contribution"
CONVENTION_B_SOURCES = ("eQuilibrator", "dGPredictor", "dGPredictor-ModelSEED")
DB_PATH_ENV = "MODELSEED_DB_PATH"

# The database marks missing estimates with 1e7
SENTINEL = 1e6
GAS_CONSTANT = 1.9858775 / 1000  # Kcal/(K mol)

# Compartment defaults of modelseedpy's FullThermoPkg: the extracellular space
# is the 0 mV reference at pH 6.5, the cytosol is at pH 7 and -160 mV, and
# concentrations lie between 1 uM and 20 mM. They are keyed by compartment
# letter, so that 'c' applies to 'c0', 'c1', ...
DEFAULT_PH = {"c": 7.0, "e": 6.5}
DEFAULT_POTENTIAL = {"c": -160.0}  # mV
DEFAULT_IONIC_STRENGTH = {"c": 0.25}  # M
DEFAULT_C_MIN = 1e-6  # M
DEFAULT_C_MAX = 0.02  # M

CPD_REGEX = re.compile(r"^(cpd\d+)")
ELEMENT_REGEX = re.compile(r"([A-Z][a-z]*)(\d*)")

logger = logging.getLogger(__name__)


def find_modelseed_database(db_path=None):
    """Locate a local clone of the ModelSEEDDatabase.

    :param str db_path: Path of the clone. Defaults to the ``MODELSEED_DB_PATH``
        environment variable, then to the ``[biochem] path`` of modelseedpy's
        configuration, read relative to the working directory and to the
        modelseedpy repository.
    :return: The absolute path of the clone
    :rtype: str
    """
    if db_path is not None:
        candidates = [db_path]
    else:
        candidates = []
        if os.environ.get(DB_PATH_ENV):
            candidates.append(os.environ[DB_PATH_ENV])
        try:
            from modelseedpy.helpers import config, project_dir

            path = config.get("biochem", "path")
            candidates += [path, os.path.join(project_dir, os.pardir, path)]
        except ImportError:
            pass

    for path in candidates:
        if os.path.isfile(os.path.join(path, "Biochemistry", "compound_00.json")):
            return os.path.abspath(path)

    raise FileNotFoundError(
        "ModelSEEDDatabase not found in {}. Pass db_path or set {}.".format(
            candidates, DB_PATH_ENV
        )
    )


@lru_cache(maxsize=2)
def _load_compounds(db_path):
    """Read the compounds of the database that have a group contribution energy.

    :return: the pyTFA metabolite entries, paired with the group contribution
        species, and the canonical (formula, charge) of the compounds, both
        keyed by ModelSEED id
    """
    biochem = os.path.join(db_path, "Biochemistry")

    species = {}
    tsv_path = os.path.join(
        biochem, "Thermodynamics", "ModelSEED", "ModelSEED_GroupContribution.tsv"
    )
    with open(tsv_path) as fid:
        rows = (
            line.rstrip("\n").split("\t") for line in fid if not line.startswith("#")
        )
        header = next(rows)
        for row in rows:
            entry = dict(zip(header, row))
            if entry["status"] == "ok":
                # The charge of the analysed species closes the groups string
                charge = int(entry["groups"].rsplit("|", 1)[-1])
                species[entry["compound_id"]] = (entry["formula"], charge)

    compounds = {}
    canonical = {}
    for filename in sorted(glob.glob(os.path.join(biochem, "compound_*.json"))):
        with open(filename) as fid:
            records = json.load(fid)
        for record in records:
            energy = (record.get("thermodynamics") or {}).get(GROUP_CONTRIBUTION)
            if not energy or energy[0] >= SENTINEL or energy[1] >= SENTINEL:
                continue
            canonical[record["id"]] = (record["formula"], record["charge"])
            # Compounds injected without groups (H+, H2O, ...) are given for
            # their canonical species
            formula, charge = species.get(record["id"], canonical[record["id"]])
            mass = record.get("mass")
            compounds[record["id"]] = {
                "id": record["id"],
                "name": record["name"],
                "formula": formula,
                "charge_std": charge,
                "nH_std": count_atoms(formula)["H"],
                "mass_std": BIGM_DG if mass is None else mass,
                "deltaGf_std": energy[0],
                "deltaGf_err": energy[1],
                "error": "Nil",
                "pKa": get_acidic_pkas(record),
                "struct_cues": {},
                "other_names": [],
            }
    return compounds, canonical


def count_atoms(formula):
    """Number of atoms of each element in a chemical formula"""
    atoms = Counter()
    for element, number in ELEMENT_REGEX.findall(formula or ""):
        atoms[element] += int(number) if number else 1
    return atoms


def get_acidic_pkas(record):
    """Acidic pKas of a ModelSEED compound record, sorted in ascending order.

    Marvin's pKas are strings like ``'1:12.60;1:3.33'``. The legacy ``pka``
    field, with an extra atom index, is used when there is no Marvin entry.
    """
    pkas = record.get("pkas") or {}
    text = pkas["Marvin"]["pKa"] if "Marvin" in pkas else record.get("pka")
    return sorted(
        float(item.rsplit(":", 1)[-1]) for item in (text or "").split(";") if item
    )


def _least_protonated_species(entry, temperature, min_ph, max_ph):
    """Energy, charge and number of protons of the least protonated species
    pyTFA starts its transform from (:meth:`MetaboliteThermo.calcDGspA`)"""
    # A MetaboliteThermo replaces its attributes by its values once initialized,
    # so the method is called on an instance holding only what it reads
    thermo = MetaboliteThermo.__new__(MetaboliteThermo)
    thermo.__dict__.update(
        debug=False,
        id=entry["id"],
        pKa=entry["pKa"],
        charge_std=entry["charge_std"],
        nH_std=entry["nH_std"],
        deltaGf_std=entry["deltaGf_std"],
        RT=GAS_CONSTANT * temperature,
        MIN_pH=min_ph,
        MAX_pH=max_ph,
    )
    return thermo.calcDGspA()


def pair_with_species(
    entry,
    formula,
    charge,
    temperature=std.TEMPERATURE_0,
    min_ph=std.MIN_PH,
    max_ph=std.MAX_PH,
):
    """Pair a metabolite entry with another protonation state of the compound.

    The energy is shifted so that pyTFA reaches the same least protonated
    species from both, which leaves the transformed energy of formation
    unchanged at any pH and ionic strength.

    :param dict entry: a thermo_data metabolite entry
    :param str formula: formula of the other species
    :param int charge: charge of the other species
    :param temperature, min_ph, max_ph: those of the ThermoModel
    :return: the new entry, or None if pyTFA cannot relate the two species
    :rtype: dict
    """
    if formula == entry["formula"] and charge == entry["charge_std"]:
        return entry
    if formula is None or charge is None:
        return None

    atoms = count_atoms(formula)
    entry_atoms = count_atoms(entry["formula"])
    protons = atoms.pop("H", 0)
    entry_protons = entry_atoms.pop("H", 0)
    if atoms != entry_atoms or protons - entry_protons != charge - entry["charge_std"]:
        return None

    new = dict(entry, formula=formula, charge_std=charge, nH_std=protons, deltaGf_std=0)
    energy, *species = _least_protonated_species(entry, temperature, min_ph, max_ph)
    offset, *new_species = _least_protonated_species(new, temperature, min_ph, max_ph)
    if species != new_species:
        return None
    new["deltaGf_std"] = energy - offset
    return new


def get_seed_id_candidates(metabolite):
    """ModelSEED ids of a metabolite, from its id then its annotations"""
    candidates = []
    values = [metabolite.id]
    for key in ["seed_id", "seed.compound"] + sorted(metabolite.annotation):
        value = metabolite.annotation.get(key)
        values += value if isinstance(value, list) else [value]
    for value in values:
        match = CPD_REGEX.match(value) if isinstance(value, str) else None
        if match and match[1] not in candidates:
            candidates.append(match[1])
    return candidates


def annotate_seed_ids(model, db_path=None, overwrite=False):
    """Annotate the metabolites with the ``seed_id`` pyTFA uses to find their
    thermodynamic data.

    When a metabolite has several ModelSEED ids (e.g. the ``seed.compound``
    annotation of BiGG models), the first one with a group contribution energy
    is chosen. The other annotations are kept.

    :param model: cobra.Model, or a modelseedpy MSModelUtil
    :param str db_path: see :func:`find_modelseed_database`
    :param bool overwrite: replace existing ``seed_id`` annotations
    :return: ids of the metabolites without any ModelSEED id
    :rtype: list(str)
    """
    model = getattr(model, "model", model)
    compounds, _ = _load_compounds(find_modelseed_database(db_path))

    missing = []
    for met in model.metabolites:
        if "seed_id" in met.annotation and not overwrite:
            continue
        candidates = get_seed_id_candidates(met)
        if not candidates:
            missing.append(met.id)
            continue
        met.annotation["seed_id"] = next(
            (x for x in candidates if x in compounds), candidates[0]
        )
    return missing


def build_thermo_from_modelseed(
    model=None,
    db_path=None,
    source=GROUP_CONTRIBUTION,
    temperature=std.TEMPERATURE_0,
    min_ph=std.MIN_PH,
    max_ph=std.MAX_PH,
    keep_unpaired=True,
):
    """Build the `thermo_data` structure from the ModelSEEDDatabase.

    The structure of the returned dictionary is specified in the pyTFA
    [documentation](https://pytfa.readthedocs.io/en/latest/thermoDB.html).

    :param model: cobra.Model or modelseedpy MSModelUtil, annotated with
        :func:`annotate_seed_ids`. If given, only its compounds are kept, and
        they are paired with the species of its metabolites.
    :param str db_path: see :func:`find_modelseed_database`
    :param str source: the ``thermodynamics`` source to read energies from
    :param temperature, min_ph, max_ph: those the ThermoModel will be built with
    :param bool keep_unpaired: keep the compounds that cannot be paired with the
        species of the model. :meth:`ThermoModel.prepare` then replaces their
        formulas and may add protons to their reactions, which changes the
        model. If False, they are left out and their reactions are not
        thermodynamically constrained.
    :return thermo_data: dict
        to be passed as argument to initialize a `ThermoModel`.
    """
    if source != GROUP_CONTRIBUTION:
        if source in CONVENTION_B_SOURCES:
            raise NotImplementedError(
                "{} energies are already transformed to pH 7 (Convention B), "
                "pyTFA would transform them again. Use '{}'.".format(
                    source, GROUP_CONTRIBUTION
                )
            )
        raise ValueError("Unknown ModelSEED thermodynamics source: " + source)

    compounds, canonical = _load_compounds(find_modelseed_database(db_path))

    # Species to pair each compound with, by order of preference
    targets = {k: [v] for k, v in canonical.items()}
    if model is not None:
        model = getattr(model, "model", model)
        model_species = {}
        for met in model.metabolites:
            if "seed_id" in met.annotation:
                model_species.setdefault(met.annotation["seed_id"], set()).add(
                    (met.formula, met.charge)
                )
        if not model_species:
            raise ValueError(
                "No metabolite has a seed_id annotation, see annotate_seed_ids"
            )
        compounds = {k: v for k, v in compounds.items() if k in model_species}
        for seed_id in compounds:
            # A compound with several species in the model is left canonical
            if len(model_species[seed_id]) == 1:
                targets[seed_id].insert(0, next(iter(model_species[seed_id])))

    metabolites = {}
    for seed_id, entry in compounds.items():
        paired = (
            pair_with_species(entry, formula, charge, temperature, min_ph, max_ph)
            for formula, charge in targets[seed_id]
        )
        metabolites[seed_id] = next((x for x in paired if x is not None), entry)

    if model is not None:
        unpaired = sorted(
            k
            for k, v in metabolites.items()
            if (v["formula"], v["charge_std"]) not in model_species[k]
        )
        if unpaired and keep_unpaired:
            logger.warning(
                "%d compounds keep a protonation state differing from the "
                "model's, and ThermoModel.prepare() may add protons to their "
                "reactions: %s",
                len(unpaired),
                ", ".join(unpaired),
            )
        elif unpaired:
            logger.warning(
                "%d compounds are left out as their protonation state differs "
                "from the model's, and their reactions get no thermodynamic "
                "constraints: %s",
                len(unpaired),
                ", ".join(unpaired),
            )
            for seed_id in unpaired:
                del metabolites[seed_id]

    return {
        "name": "ModelSEED " + source,
        "units": "kcal/mol",
        "cues": {},
        "metabolites": deepcopy(metabolites),
    }


def _lookup(value, symbol, defaults, fallback):
    """Value of a compartment parameter given as a number, or as a dict keyed by
    compartment symbol or letter"""
    letter = symbol.rstrip("0123456789")
    if value is None:
        value = {}
    if not isinstance(value, dict):
        return value
    for key in (symbol, letter):
        if key in value:
            return value[key]
    return defaults.get(letter, fallback)


def build_compartment_data(
    model,
    ph=None,
    ionic_strength=None,
    potential=None,
    c_min=DEFAULT_C_MIN,
    c_max=DEFAULT_C_MAX,
):
    """Build the compartment data of a model, to be applied with
    :func:`pytfa.io.apply_compartment_data`.

    Each parameter is a number applying to every compartment, or a dict keyed by
    compartment symbol (``'c0'``) or letter (``'c'``). Unlisted compartments get
    the defaults of modelseedpy's FullThermoPkg, and else pH 7, 0 M ionic
    strength and 0 mV.

    :param model: cobra.Model or modelseedpy MSModelUtil
    :param ph: pH
    :param ionic_strength: ionic strength (M)
    :param potential: electrical potential (mV) relative to the extracellular
        space
    :param c_min: minimal metabolite concentration (M)
    :param c_max: maximal metabolite concentration (M)
    :return: dict
    """
    model = getattr(model, "model", model)
    symbols = sorted({met.compartment for met in model.metabolites})
    potentials = {x: _lookup(potential, x, DEFAULT_POTENTIAL, 0.0) for x in symbols}

    compartment_data = {}
    for symbol in symbols:
        name = model.compartments.get(symbol)
        if isinstance(name, dict):
            name = name.get("name")
        compartment_data[symbol] = {
            "symbol": symbol,
            "name": name or symbol,
            "pH": _lookup(ph, symbol, DEFAULT_PH, 7.0),
            "ionicStr": _lookup(ionic_strength, symbol, DEFAULT_IONIC_STRENGTH, 0.0),
            "c_min": _lookup(c_min, symbol, {}, DEFAULT_C_MIN),
            "c_max": _lookup(c_max, symbol, {}, DEFAULT_C_MAX),
            # Potential of each compartment as seen from this one, as used for
            # transport from this compartment to the other
            "membranePot": {x: potentials[x] - potentials[symbol] for x in symbols},
        }
    return compartment_data
