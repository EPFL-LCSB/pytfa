# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Thermodynamic computations for reactions


"""
from functools import reduce
from math import log, sqrt

from . import std
from .utils import find_transported_mets
from .metabolite import CPD_PROTON

###################
# REACTIONS TOOLS #
###################



def calcDGtpt_rhs(reaction, compartmentsData, thermo_units):
    """ Calculates the RHS of the deltaG constraint, i.e. the sum of the
    non-concentration terms

    :param cobra.thermo.reaction.Reaction reaction: The reaction to compute the
        data for
    :param dict(float) compartmentsData: Data of the compartments of the cobra_model
    :param str thermo_units: The thermodynamic database of the cobra_model

    :returns: deltaG_tpt and the breakdown of deltaG_tpt
    :rtype: tuple(float, dict(float))

    Example:
        ATP Synthase reaction::

            reaction = cpd00008 + 4 cpd00067 + cpd00009 <=> cpd00002 + 3 cpd00067 + cpd00001
            compartments =  'c'       'e'        'c'           'c'         'c'         'c'

    If there are any metabolites with unknown energies then returns
        ``(0, None)``.

    """

    # Compute our constants in accordance with the thermoDB
    if thermo_units == "kJ/mol":
        GAS_CONSTANT = 8.314472 / 1000  # kJ/(K mol)
        faraday_const = 96.485  # kJ/eV
    else:
        GAS_CONSTANT = 1.9858775 / 1000  # Kcal/(K mol)
        faraday_const = 23.061  # kcal/eV

    TEMPERATURE = 298.15  # K

    RT = GAS_CONSTANT * TEMPERATURE

    if reduce(lambda count, met: (
                count + (1 if met.thermo.deltaGf_tr > 10 ** 6 else count)),
              reaction.metabolites,
              0) > 1:
        return (0, None)

    sum_deltaGFis_trans = 0
    sum_stoich_NH = 0
    RT_sum_H_LC_tpt = 0  # to include the differential proton concentration
    # effects if protons are transported

    transportedMets = find_transported_mets(reaction)
    compartments = {'reactant': [], 'product': []}

    for seed_id in transportedMets:
        for metType in ['reactant', 'product']:
            if seed_id != 'cpd00001':
                met = transportedMets[seed_id][metType]
                pH_comp = met.thermo.pH
                ionicStr_comp = met.thermo.ionicStr

                deltaGfsp = met.thermo.deltaGf_tr

                compartments[metType].append(met.compartment)
                sum_stoich_NH += ((1 if metType == 'product' else -1)
                                  * transportedMets[seed_id]['coeff']
                                  * met.thermo.nH_std
                                  * RT
                                  * log(10 ** -pH_comp))
                sum_deltaGFis_trans += ((1 if metType == 'product' else -1)
                                        * transportedMets[seed_id]['coeff']
                                        * deltaGfsp)
            else:
                compartments[metType].append('')

            if seed_id == CPD_PROTON:
                met = transportedMets[seed_id][metType]
                pH_comp = met.thermo.pH
                RT_sum_H_LC_tpt += ((1 if metType == 'product' else -1)
                                    * RT
                                    * transportedMets[seed_id]['coeff']
                                    * log(10 ** -pH_comp))

    # calculate the transport of any ions
    # membrane potential is always defined as inside - outside
    # we should take the larger stoich of the transported compound
    sum_F_memP_charge = 0

    for seed_id in transportedMets:
        if seed_id != 'cpd00001':
            out_comp = transportedMets[seed_id]['reactant'].compartment
            in_comp = transportedMets[seed_id]['product'].compartment
            mem_pot = compartmentsData[out_comp]['membranePot'][in_comp]
            charge = transportedMets[seed_id]['reactant'].thermo.charge_std
            # Equal to the product's one
            sum_F_memP_charge += (faraday_const
                                  * (mem_pot / 1000.)
                                  * transportedMets[seed_id]['coeff']
                                  * charge)

    deltaG = 0

    for met in reaction.metabolites:
        if CPD_PROTON != met.annotation['seed_id']:
            deltaG += reaction.metabolites[met] * met.thermo.deltaGf_tr

    sum_deltaGFis = 0

    # lastly we calculate the deltaG of the chemical reaction if any
    # but we do not add this part to the rhs as it would be included in the
    # potential energy of the enzyme


    final_coeffs = reaction.metabolites.copy()

    for seed_id in transportedMets:
        for metType in ['reactant', 'product']:
            final_coeffs[transportedMets[seed_id][metType]] -= (
                (1 if metType == 'product' else -1)
                * transportedMets[seed_id]['coeff'])

    for met in final_coeffs:
        if final_coeffs[met] != 0 and met.annotation['seed_id'] != CPD_PROTON:

            met_deltaGis = met.thermo.deltaGf_tr
            sum_deltaGFis += final_coeffs[met] * met_deltaGis

    # Sum all the parts
    DG_trans_RHS = (sum_stoich_NH
                    + sum_F_memP_charge
                    + sum_deltaGFis_trans
                    + RT_sum_H_LC_tpt
                    + sum_deltaGFis)

    breakdown = {
        'sum_deltaGFis': sum_deltaGFis,
        'sum_stoich_NH': sum_stoich_NH,
        'sum_F_memP_charge': sum_F_memP_charge,
        'sum_deltaGFis_trans': sum_deltaGFis_trans,
        'RT_sum_H_LC_tpt': RT_sum_H_LC_tpt
    }

    return (DG_trans_RHS, breakdown)


def calcDGR_cues(reaction, reaction_cues_data):
    """ Calculates the deltaG reaction and error of the reaction using the
    constituent structural cues changes and returns also the error if any.

    :param cobra.thermo.reaction.Reaction reaction: The reaction to compute
        deltaG for
    :param dict reaction_cues_data:

    :returns: deltaGR, error on deltaGR, the cues in the reaction (keys of the
        dictionnary) and their indices (values of the dictionnary),
        and the error code if any.

        If everything went right, the error code is an empty string

    :rtype: tuple(float, float, dict(float), str)

    """

    deltaGR = 0
    deltaGR_err = 0
    cues = {}
    error = ''

    # First we should check if all the reactants are in terms of compound IDs
    for reactant in reaction.metabolites:
        if len(reactant.thermo.struct_cues) == 0:
            return (10 ** 7, 10 ** 7, '', 'UNKNOWN_GROUPS')
        (deltaGF, deltaGFerr, cpd_cues) = calcDGF_cues(
            reactant.thermo.struct_cues,
            reaction_cues_data)
        for cue in cpd_cues:
            if cue in cues:
                cues[cue] += reaction.metabolites[reactant] * cpd_cues[cue]
            else:
                cues[cue] = reaction.metabolites[reactant] * cpd_cues[cue]

    for cue in cues:
        deltaGR += cues[cue] * reaction_cues_data[cue]['energy']
        deltaGR_err += (cues[cue] * reaction_cues_data[cue]['error']) ** 2

    deltaGR_err = sqrt(deltaGR_err)

    return (deltaGR, deltaGR_err, cues, error)


def calcDGF_cues(cues, reaction_cues_data):
    """ Calculates the deltaG formation and error of the compound using its
    constituent structural cues.

    :param list(str) cues: A list of cues' names
    :param dict reaction_cues_data:

    :returns: deltaG formation, the error on deltaG formation, and a dictionnary
        with the cues' names as key and their coefficient as value
    :rtype: tuple(float, float, dict(float)).

    """
    deltaGF = 0
    deltaGF_err = 0
    finalcues = {}

    for cue in cues:
        if cue in finalcues:
            finalcues[cue] += cues[cue]
        else:
            finalcues[cue] = cues[cue]

        deltaGF += reaction_cues_data[cue]['energy'] * cues[cue]
        deltaGF_err += (reaction_cues_data[cue]['error'] * cues[cue]) ** 2

    deltaGF_err = sqrt(deltaGF_err)

    return (deltaGF, deltaGF_err, finalcues)


def get_debye_huckel_b(T):
    """
    The Debye-Huckel A and B do depend on the temperature
    As for now though they are returned as a constant (value at 298.15K)

    :param T: Temperature in Kelvin
    :return: Debye_Huckel_B
    """
    return std.DEBYE_HUCKEL_B_0