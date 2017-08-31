# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Thermodynamic computations


"""

from math import log, sqrt
from functools import reduce

from .utils import find_transported_mets
from . import std

from ..utils.numerics import BIGM_DG

CPD_PROTON = 'cpd00067'

DEFAULT_VAL = BIGM_DG

class MetaboliteThermo:
    """
    A class representing the thermodynamic values of a metabolite

    :param dict metData: A dictionnary containing the values for the
        metabolite, from the thermodynamic database
    :param float pH: The pH of the metabolite's compartment
    :param float ionicStr: The ionic strength of the metabolite's
        compartment
    :param temperature:
    :param min_ph:
    :param max_ph:
    :param debye_huckel_b:
    :param string thermo_unit: The unit used in `metData`'s values
    :param bool debug: *Optional* If set to ``True``, some debugging values
        will be printed. This is only useful for development or debugging
        purposes.

    .. note::

        The values are automatically computed on class creation. Usually you
        don't have to call any methods defined by this class, but only to access
        the attributes you need.

    The available attributes are :

    +-------------+-----------------------------------------------------------+
    | id          | The metabolite's seed_id                                   |
    +-------------+-----------------------------------------------------------+
    | pKa         | pKas of the metabolite                                    |
    +-------------+-----------------------------------------------------------+
    | error       | Error on the metabolite's thermodynamic data. This is the |
    |             | value from the thermodynamic database.                    |
    |             | Thermodynamic values will be computed only if this equals |
    |             | to 'Nil'.                                                 |
    +-------------+-----------------------------------------------------------+
    | deltaGf_std | Transformed Gibbs energy of formation of the metabolite,  |
    |             | in standard conditions (value from the the thermodynamic  |
    |             | database).                                                |
    +-------------+-----------------------------------------------------------+
    | deltaGf_err | Error on the transformed Gibbs energy of formation of the |
    |             | metabolite,in standard conditions (value from the the     |
    |             | thermodynamic database).                                  |
    +-------------+-----------------------------------------------------------+
    | mass        | Mass of the metabolite (value from the thermodynamic      |
    |             | database).                                                |
    +-------------+-----------------------------------------------------------+
    | nH_std      | Number of protons of the metabolite, in standard          |
    |             | conditions (extracted from the thermodynamic database).   |
    +-------------+-----------------------------------------------------------+
    | struct_cues | The cues that make the structure of the metabolite.       |
    +-------------+-----------------------------------------------------------+
    | charge_std  | The charge (mV) of the metabolite, in standard conditions.|
    |             | This value is extracted from the thermodynamic database.  |
    +-------------+-----------------------------------------------------------+
    | pH          | The pH of the metabolite's compartment                    |
    +-------------+-----------------------------------------------------------+
    | ionicStr    | The ionic strength of the metabolite's compartment        |
    +-------------+-----------------------------------------------------------+
    | deltaGf_tr  |  Transformed Gibbs energy of formation of specie with     |
    |             |  given pH and ionic strength using formula given by       |
    |             |  Goldberg and Tewari, 1991                                |
    +-------------+-----------------------------------------------------------+

    Since the reactions expose similar values through a dictionnary, it is
    better to access the attributes aforementionned of this class as if it was
    a dictionnary : ``metabolite.thermo['pH']``.

    """

    def __init__(self, metData, pH, ionicStr, temperature=std.TEMPERATURE_0,
                 min_ph=std.MIN_PH, max_ph=std.MAX_PH,
                 debye_huckel_b=std.DEBYE_HUCKEL_B_0, thermo_unit='kJ/mol',
                 debug=False):
        # Store all the data we got
        """

        :param metData:
        :param pH:
        :param ionicStr:
        :param temperature:
        :param min_ph:
        :param max_ph:
        :param debye_huckel_b:
        :param thermo_unit:
        :param debug:
        """
        self.debug = debug
        if self.debug and metData != None:
            print("Debugging thermo computations for metabolite : {} ({})"
                  .format(metData['name'], metData['id'])
                  )

        self.pH = pH
        self.ionicStr = ionicStr
        self.Debye_Huckel_B = debye_huckel_b

        # Compute internal values to adapt the the thermo_unit provided
        if thermo_unit == "kJ/mol":
            GAS_CONSTANT = 8.314472 / 1000  # kJ/(K mol)
            self.Adjustment = 1
        else:
            GAS_CONSTANT = 1.9858775 / 1000  # Kcal/(K mol)
            self.Adjustment = 4.184

        TEMPERATURE = temperature

        self.RT = GAS_CONSTANT * TEMPERATURE

        # CONSTANTS
        self.MAX_pH = max_ph
        self.MIN_pH = min_ph

        # Store metData values
        self.id = None if metData == None else metData['id']
        self.pKa = [] if metData == None else metData['pKa']
        self.error = None if metData == None else metData['error']
        self.deltaGf_std = DEFAULT_VAL if metData == None else metData[
            'deltaGf_std']
        self.deltaGf_err = DEFAULT_VAL if metData == None else metData[
            'deltaGf_err']
        self.mass = DEFAULT_VAL if metData == None else metData['mass_std']
        self.nH_std = None if metData == None else metData['nH_std']
        self.struct_cues = [] if metData == None else metData['struct_cues']
        self.charge_std = DEFAULT_VAL if metData == None else metData['charge_std']
        self.struct_cues = None if metData == None else metData['struct_cues']

        # Compute deltaGf_tr if possible
        self.deltaGf_tr = DEFAULT_VAL if metData == None else self.calcDGis()

        self.__dict__ = {
            'id': self.id,
            'pKa': self.pKa,
            'error': self.error,
            'deltaGf_std': self.deltaGf_std,
            'deltaGf_err': self.deltaGf_err,
            'mass': self.mass,
            'nH_std': self.nH_std,
            'charge_std': self.charge_std,
            'struct_cues': self.struct_cues,
            'deltaGf_tr': self.deltaGf_tr,
            'pH': self.pH,
            'ionicStr': self.ionicStr,
        }

    # Various methods to have a dictionnary-like behavior, for consistency with
    # the reactions' thermo attribute
    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self):
        return repr(self.__dict__)

    def keys(self):
        return self.__dict__.keys()

    def values(self):
        return self.__dict__.values()

    def items(self):
        return self.__dict__.items()

    def __cmp__(self, dict_):
        return cmp(self.__dict__, dict_)

    def __contains__(self, item):
        return item in self.__dict__

    def __iter__(self):
        return iter(self.__dict__)

    def __unicode__(self):
        return unicode(repr(self.__dict__))

    def calcDGis(self):
        """ Calculate the transformed Gibbs energy of formation of specie with
        given pH and ionic strength using formula given by Goldberg and Tewari,
        1991

        Equation 4.5-6 in Alberty's book

        :returns: DG_is for the metabolite
        :rtype: float

        """

        if self.debug:
            print("Computing DGis...")

        # Special case for protons...
        if self.id == CPD_PROTON:
            if self.debug:
                print("Found proton")
            return -self.RT * log(10 ** -self.pH)

        # Error too big
        if self.error != 'Nil' or self.deltaGf_std > 9 * 10 ** 6:
            if self.debug:
                print("Error or deltaGf_std too big")
            return 10 ** 7

        # Compute the value
        P = self.calc_potential()
        DGis = self.calcDGsp() - self.RT * log(P)

        if self.debug:
            print("Found DGis = " + str(DGis))

        return DGis

    def calcDGsp(self):
        """ Calculate the transformed Gibbs energy of formation of specie with
        given pH and ionic strength using formula given by Goldberg and Tewari,
        1991

        Equation 4.4-10 in Alberty's book

        :returns: DG_sp for the metabolite
        :rtype: float


        """
        (deltaGo, charge, nH) = self.calcDGspA()
        zsq = charge ** 2
        I = self.ionicStr
        term1 = (nH * self.RT * log(10 ** -self.pH))
        term2 = (2.91482
                 * (zsq - nH)
                 * sqrt(I)
                 / (1 + self.Debye_Huckel_B * sqrt(I))
                 ) / self.Adjustment

        return deltaGo - (term1 + term2)

    def calc_potential(self):
        """ Calculate the binding polynomial of a specie, with pK values

        :returns: The potential of the metabolite
        :rtype: float

        """
        if self.debug:
            print("Computing P...")

        # Init some values used here...
        prod_denom = 1
        p = 1
        pka_values = self.get_pka()

        if self.debug:
            print("pKas used : " + str(pka_values))

        # Make the computation...
        if len(pka_values) > 0:
            if min(pka_values) <= self.MAX_pH:
                for i, this_pka in enumerate(pka_values):
                    numerator = 10 ** (-(i + 1) * self.pH)
                    K = 10 ** (-this_pka)
                    denominator = prod_denom * K
                    prod_denom = denominator
                    if self.debug:
                        print("Adding "
                              + str(numerator / denominator)
                              + " to P")

                    p += numerator / denominator

        return p

    def get_pka(self):
        """ Get the pKas of the metabolite

        :returns: The pKas of the metabolite
        :rtype: list(float)

        """
        if self.debug:
            print("Getting the list of pKas...")
        (deltaGspA, charge, sp_nH) = self.calcDGspA()

        pka_list = self.pKa

        pka_values = [None] * len(pka_list)

        # Get only useful pKas
        pka_list = [x for x in pka_list if 3 < x < 9]
        # Sort the list
        pka_list.sort(reverse=True)

        j = 0

        if len(pka_list) > 1:
            for i,this_pka in enumerate(pka_list):
                sigmanusq = 1 + (charge + i) ** 2 - (charge + i - 1) ** 2
                if self.MAX_pH > pka_list[i] > self.MIN_pH:
                    pka_values[j] = self._calc_pka(this_pka,sigmanusq)

                    if self.debug:
                        print("Added to pKas : " + str(pka_values[j]))

                    j += 1

        elif len(pka_list) == 1:
            if self.debug:
                print("Only one pKa to add")

            sigmanusq = 2 * charge
            pka_values[j] = self._calc_pka(pka_list[j],sigmanusq)

        # Only return useful values
        pka_values = [x for x in pka_values if x != None]

        if self.debug:
            print("Filtered pKa values : " + str(pka_values))

        pka_values.sort(reverse=True)

        return pka_values

    def _calc_pka(self, pka,sigmanusq):
        lnkzero = log(10 ** -pka)
        pka_value = -(
            lnkzero - sigmanusq * (1.17582 * sqrt(self.ionicStr)) / (
            1 + 1.6 * sqrt(self.ionicStr))) / log(10)
        return pka_value

    def calcDGspA(self):
        """ Calculates deltaGf, charge and nH of the specie when it is at least
        protonated state based on MFAToolkit compound data for the pKa values
        within the range considered (MIN_pH to MAX_pH).

        These values are used as the starting point for Alberty's calculations.

        :returns: deltaGspA, sp_charge and sp_nH
        :rtype: tuple(float, float, int)

        """
        if self.debug:
            print('Computing DGspA()...')

        # Case of the proton
        if self.id == CPD_PROTON:
            if self.debug:
                print('Proton found, returning standard values')
            # we do not adjust for proton so just return the values
            deltaGspA = self.deltaGf_std
            sp_charge = self.charge_std
            sp_nH = self.nH_std
            return (deltaGspA, sp_charge, sp_nH)

        pka_list = self.pKa

        acceptedpKas = [x for x in pka_list if (
            self.MIN_pH < x < self.MAX_pH
        )]

        # No pKas found
        if len(acceptedpKas) == 0:
            if self.debug:
                print('No pKas found, returning standard values')
            deltaGspA = self.deltaGf_std
            sp_charge = self.charge_std
            sp_nH = self.nH_std
            return (deltaGspA, sp_charge, sp_nH)

        # if compound is multiprotic
        # discard pKa values above MAX_pH

        pka_list = [x for x in pka_list if x < self.MAX_pH]
        charge_adj = reduce(lambda count, x: count + 1 if x < 7 else count,
                            acceptedpKas,
                            0)
        sp_charge = -len(pka_list)

        pKs = []
        sp_deltaGf = 0
        sp_nH = 0

        # if the ReactionDB compound is already at the right state, i.e. most
        # negative state we can just return their values too after adjusting
        # for the ionic strength

        if self.charge_std == sp_charge:
            if self.debug:
                print('charge found is equal to standard charge, '
                      + 'returning default values')
            deltaGspA = self.deltaGf_std
            sp_charge = self.charge_std
            sp_nH = self.nH_std
            return (deltaGspA, sp_charge, sp_nH)

        charge = self.charge_std

        num_iter = charge - sp_charge

        sp_nH = self.nH_std - num_iter
        deltaGspA = self.deltaGf_std

        npKa = len(pka_list)

        if npKa - num_iter == -1:
            start = 0  # NOTE : In Python, lists begin at 0
        else:
            start = npKa - num_iter

        pKs = pka_list[start:]

        if self.debug:
            print("Accepted pKs : ", pKs)

        if self.debug:
            print(pka_list, start, pKs)

        for j,_ in enumerate(pKs):
            deltaGspA -= self.RT * log(10 ** -pKs[j])

        if self.debug:
            print("Found deltaGspA : " + str(deltaGspA))

        return (deltaGspA, sp_charge, sp_nH)


###################
# REACTIONS TOOLS #
###################

def calcDGtpt_rhs(reaction, compartmentsData, thermo_units):
    """ Calculates the RHS of the deltaG constraint, i.e. the sum of the
    non-concentration terms

    :param cobra.core.reaction.Reaction reaction: The reaction to compute the
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
    # potential energy of the metabolite


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

    :param cobra.core.reaction.Reaction reaction: The reaction to compute
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