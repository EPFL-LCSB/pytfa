"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Some tools around COBRApy models used by pyTFA


"""
import re

Formula_regex = re.compile("([A-Z][a-z]*)([0-9]*)")

# idenfifiers of water metabolite in different databases
PROTON = {
    "h",  "15378", "10744", "13357", "5584", "24636", "29233", "29234",
    "MNXM104313", "MNXM113751", "MNXM145872", "MNXM89553", "HMDB59597",
    "C00080", "PROTON", "1132304", "113529", "1470067", "156540", "163953",
    "193465", "194688", "2000349", "2872447", "351626", "372511", "374900",
    "425969", "425978", "425999", "427899", "428040", "428548", "5668577",
    "70106", "74722",  "39",  "cpd00067",  "MNXM1"
}

# idenfifiers of water metabolite in different databases
WATER = {
    "h2o", "oh1", "15377", "10743", "13352", "27313", "42043", "42857",
    "43228", "44292", "44701", "44819", "5585", "16234", "13365", "13419",
    "44641", "5594", "29356", "29374", "29375", "29412", "30490", "33806",
    "33811", "33813", "41981", "29373", "41979", "MNXM114710", "MNXM114753",
    "MNXM11838", "MNXM124004", "MNXM124324", "MNXM124831", "MNXM125045",
    "MNXM126600", "MNXM128935", "MNXM131091", "MNXM145357", "MNXM49218",
    "MNXM56889", "MNXM89551", "5882df9c-dae1-4d80-a40e-db4724271456\/compound\/969d0227-3069-4e44-9525-7ae7bad84170",
    "650babc9-9d68-4b73-9332-11972ca26f7b\/compound\/799908db-b8c9-4982-86cb-1f225e2ad08c",
    "650babc9-9d68-4b73-9332-11972ca26f7b\/compound\/e7f34a8e-cded-4793-b6d5-792335b38636",
    "HMDB02111", "C00001", "D00001", "C01328", "C18714", "D03703", "D06322",
    "CPD-15815", "HYDROXYL-GROUP", "OH", "OXONIUM", "WATER", "109276",
    "1130930", "113518", "113519", "113521", "141343", "1605715", "189422",
    "2022884", "29356", "351603", "5278291", "5668574", "5693747", "8851517",
    "40", "cpd00001", "cpd15275", "cpd27222", "MNXM2"
}


def check_reaction_balance(reaction, proton = None):
    """
    Check the balance of a reaction, and eventually add protons to balance
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
                # Jorge Carrasco: Can we use warnings here?
                reaction.model.logger.warning(
                    'Warning : ' + metabolite.formula + '/' + atom[0]
                )
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


def find_transported_mets(reaction, annotation_key="seed_id"):
    """
    Get a list of the transported metabolites of the reaction.

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
        reactants_coeffs[met.annotation[annotation_key]] = {
            'coeff':reaction.metabolites[met],
            'met':met
        }

    # Initialize the dictionnary we will fill with the transported metabolites
    trans_coeffs = {}

    for met in reaction.products:
        seed_id = met.annotation[annotation_key]

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


def check_transport_reaction(reaction, annotation_key="seed_id"):
    """

    Check if a reaction is a transport reaction

    :param cobra.thermo.reaction.Reaction reaction: The reaction to check

    :returns: Whether the reaction is a transport reaction or not
    :rtype: bool

    A transport reaction is defined as a reaction that has the same compound
    as a reactant and a product. We can distinguish them thanks to their seed_ids.
    If they have one
    If not, use met_ids and check if they are the same, minus compartment

    """
    seed_ids = {}
    try:
        for reactant in reaction.reactants:
            seed_ids[reactant.annotation[annotation_key]] = True

        for product in reaction.products:
            if product.annotation[annotation_key] in seed_ids:
                return True
    except KeyError:
        reactants_ids = [x.id.replace(x.compartment,'') for x in reaction.reactants]
        product_ids = [x.id.replace(x.compartment,'') for x in reaction.products]

        return set(reactants_ids) == set(product_ids)

    return False


def is_same_stoichiometry(this_reaction, that_reaction):
    this_met_dict = {k.id:v for k,v in this_reaction.metabolites.items()}
    that_met_dict = {k.id:v for k,v in that_reaction.metabolites.items()}
    return this_met_dict == that_met_dict


def is_exchange(rxn):
    return len(rxn.metabolites) == 1


def get_reaction_compartment(reaction):
    """Get the compartment of a reaction to then prepare it for conversion."""
    comp = None
    for met in reaction.metabolites:
        if comp is None:
            comp = met.compartment
        elif met.compartment != comp:
            comp = 'c'
    return comp
