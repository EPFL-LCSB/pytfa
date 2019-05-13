"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Some tools around COBRApy models used by pyTFA


"""
import re

Formula_regex = re.compile("([A-Z][a-z]*)([0-9]*)")


def check_reaction_balance(reaction, proton = None):
    """ Check the balance of a reaction, and eventually add protons to balance
    it

    :param cobra.thermo.reaction.Reaction reaction: The reaction to check the
        balance of.
    :param cobra.thermo.metabolite.Metabolite proton: *Optional* The proton to add
        to the reaction to balance it.

    :returns: The balance of the reaction :

        * ``drain flux``
        * ``missing structures``
        * ``balanced``
        * ``N protons added to reactants`` with ``N`` a :any:`float`
        * ``N protons added to products`` with ``N`` a :any:`float`
        * ``missing atoms``

    :rtype: str

    If ``proton`` is provided, this function will try to balance the equation
    with it, and return the result.

    If no ``proton`` is provided, this function will not try to balance the
    equation.

    .. warning::
       This function does not verify if ``proton`` is in the correct
       compartment, so make sure you provide the ``proton`` belonging to the
       correct compartment !

    """
    # Sum of the total charges of both sides (reactants counted negatively,
    # products counted positively)
    sum_charge = 0

    if len(reaction.metabolites) == 1:
        return 'drain flux'

    # Count the atoms on each side
    Atoms_sum = [0] * 26

    for metabolite in reaction.metabolites:
        if metabolite.formula == 'NA' or \
           metabolite.formula is None:
            return 'missing structures'
        metCharge = metabolite.thermo.charge_std
        metCoeff = reaction.metabolites[metabolite]

        sum_charge += metCharge * metCoeff

        # Parse the formula to extract atoms
        atoms = Formula_regex.findall(metabolite.formula)

        for atom in atoms:
            try:
                id_ = ['C', 'N', 'O', 'H', 'P', 'Na', 'Mg', 'S', 'Cl', 'K',
                'Ca', 'Mn', 'Fe', 'Ni', 'Co', 'Cu', 'Zn', 'As', 'Se', 'Ag',
                'Cd', 'W', 'Hg', 'R', 'Mo', 'X'].index(atom[0])
            except ValueError:
                print('Warning : ' + metabolite.formula + '/' + atom[0])
                continue

            Atoms_sum[id_] += metCoeff * (int(atom[1]) if atom[1] else 1)

    # Find unbalanced atoms
    nonNull = [i for i in range(26) if Atoms_sum[i] != 0]

    # If everything is balanced
    if len(nonNull) == 0 and sum_charge == 0:
        return 'balanced'

    # Add protons to compensate if possible
    if (proton and len(nonNull) == 1
        and nonNull[0] == 3
        and Atoms_sum[3] == sum_charge):

        reaction.add_metabolites({proton: -sum_charge})

        if sum_charge > 0:
            return str(sum_charge) + ' protons added to reactants'
        else:
            return str(-sum_charge) + ' protons added to products'

    return 'missing atoms'


def find_transported_mets(reaction):
    """ Get a list of the transported metabolites of the reaction.

    :param cobra.thermo.reaction.Reaction reaction: The reaction to get the
        transported metabolites of

    :returns: A dictionnary of the transported metabolites.
        The index corresponds to the seed_id of the transported enzyme

        The value is a dictionnairy with the following values:

            * coeff (:any:`float`):
                The stoechiomectric coefficient of the enzyme
            * reactant (:any:`cobra.thermo.enzyme.Metabolite`):
                The reactant of the reaction corresponding to the transported
                enzyme
            * product (:any:`cobra.thermo.enzyme.Metabolite`):
                The product of the reaction corresponding to the transported
                enzyme

    A transported enzyme is defined as a enzyme which is a product and a
    reactant of a reaction. We can distinguish them thanks to their seed_ids.

    """

    # First, we get all the reactants of the reaction
    reactants_coeffs = {}

    for met in reaction.reactants:
        reactants_coeffs[met.annotation['seed_id']] = {
            'coeff':reaction.metabolites[met],
            'met':met
        }

    # Initialize the dictionnary we will fill with the transported metabolites
    trans_coeffs = {}

    for met in reaction.products:
        seed_id = met.annotation['seed_id']

        # If the seed_id also corresponds to a reactant, we add it to the result
        if seed_id in reactants_coeffs:
            trans_coeffs[seed_id] = {
                'coeff'    : max(
                                abs(reactants_coeffs[seed_id]['coeff']),
                                abs(reaction.metabolites[met])
                            ),
                'reactant' : reactants_coeffs[seed_id]['met'],
                'product'  : met
            }

    return trans_coeffs


def check_transport_reaction(reaction):
    """ Check if a reaction is a transport reaction

    :param cobra.thermo.reaction.Reaction reaction: The reaction to check

    :returns: Whether the reaction is a transport reaction or not
    :rtype: bool

    A transport reaction is defined as a reaction that has the same compound
    as a reactant and a product. We can distinguish them thanks to their seed_ids.
    If they have one

    """
    seed_ids = {}
    try:
        for reactant in reaction.reactants:
            seed_ids[reactant.annotation['seed_id']] = True

        for product in reaction.products:
            if product.annotation['seed_id'] in seed_ids:
                return True
    except KeyError:
        return None

    return False


def is_same_stoichiometry(this_reaction, that_reaction):
    this_met_dict = {k.id:v for k,v in this_reaction.metabolites.items()}
    that_met_dict = {k.id:v for k,v in that_reaction.metabolites.items()}
    return this_met_dict == that_met_dict


def is_exchange(rxn):
    return len(rxn.metabolites) == 1