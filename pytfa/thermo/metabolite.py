# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Thermodynamic computations for metabolites


"""

from functools import reduce
from math import log, sqrt

from . import std
from ..utils.numerics import BIGM_THERMO

CPD_PROTON = 'cpd00067'

DEFAULT_VAL = BIGM_THERMO


class MetaboliteThermo:
    """
    A class representing the thermodynamic values of a enzyme

    :param dict metData: A dictionary containing the values for the
        enzyme, from the thermodynamic database
    :param float pH: The pH of the enzyme's compartment
    :param float ionicStr: The ionic strength of the enzyme's
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
    | id          | The enzyme's seed_id                                   |
    +-------------+-----------------------------------------------------------+
    | pKa         | pKas of the enzyme                                    |
    +-------------+-----------------------------------------------------------+
    | error       | Error on the enzyme's thermodynamic data. This is the |
    |             | value from the thermodynamic database.                    |
    |             | Thermodynamic values will be computed only if this equals |
    |             | to 'Nil'.                                                 |
    +-------------+-----------------------------------------------------------+
    | deltaGf_std | Transformed Gibbs energy of formation of the enzyme,  |
    |             | in standard conditions (value from the the thermodynamic  |
    |             | database).                                                |
    +-------------+-----------------------------------------------------------+
    | deltaGf_err | Error on the transformed Gibbs energy of formation of the |
    |             | enzyme,in standard conditions (value from the the     |
    |             | thermodynamic database).                                  |
    +-------------+-----------------------------------------------------------+
    | mass        | Mass of the enzyme (value from the thermodynamic      |
    |             | database).                                                |
    +-------------+-----------------------------------------------------------+
    | nH_std      | Number of protons of the enzyme, in standard          |
    |             | conditions (extracted from the thermodynamic database).   |
    +-------------+-----------------------------------------------------------+
    | struct_cues | The cues that make the structure of the enzyme.       |
    +-------------+-----------------------------------------------------------+
    | charge_std  | The charge (mV) of the enzyme, in standard conditions.|
    |             | This value is extracted from the thermodynamic database.  |
    +-------------+-----------------------------------------------------------+
    | pH          | The pH of the enzyme's compartment                    |
    +-------------+-----------------------------------------------------------+
    | ionicStr    | The ionic strength of the enzyme's compartment        |
    +-------------+-----------------------------------------------------------+
    | deltaGf_tr  |  Transformed Gibbs energy of formation of specie with     |
    |             |  given pH and ionic strength using formula given by       |
    |             |  Goldberg and Tewari, 1991                                |
    +-------------+-----------------------------------------------------------+

    Since the reactions expose similar values through a dictionnary, it is
    better to access the attributes aforementionned of this class as if it was
    a dictionnary : ``enzyme.thermo['pH']``.

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
        if self.debug and metData is not None:
            print("Debugging thermo computations for enzyme : {} ({})"
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

        :returns: DG_is for the enzyme
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

        :returns: DG_sp for the enzyme
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

        :returns: The potential of the enzyme
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
        """ Get the pKas of the enzyme

        :returns: The pKas of the enzyme
        :rtype: list(float)

        """
        if self.debug:
            print("Getting the list of pKas...")
        (deltaGspA, charge, sp_nH) = self.calcDGspA()

        pka_list = self.pKa

        pka_values = [None] * len(pka_list)

        # Get only useful pKas
        pka_list = [x for x in pka_list if self.MIN_pH < x < self.MAX_pH]
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
        pka_values = [x for x in pka_values if x is not None]

        if self.debug:
            print("Filtered pKa values : " + str(pka_values))

        pka_values.sort(reverse=True)

        return pka_values

    def _calc_pka(self, pka,sigmanusq):
        lnkzero = log(10 ** -pka)
        pka_value = -(
            lnkzero - sigmanusq * (std.DEBYE_HUCKEL_A * log(10) * sqrt(self.ionicStr)) / (
            1 + self.Debye_Huckel_B * sqrt(self.ionicStr))) / log(10)
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