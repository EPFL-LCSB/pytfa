"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamic constraints for Flux-Based Analysis of reactions

.. moduleauthor:: pyTFA team

Input/Output tools to import or export pytfa models


"""

import pickle
import re  # Regular expressions
import zlib

from cobra import Model, Reaction, Metabolite
from scipy.io import loadmat


def import_matlab_model(path, variable_name = None):
    """ Convert at matlab model to a pyTFA model, with Thermodynamic values

    :param string path: The path of the file to import

    :returns: The converted model
    :rtype: cobra.core.model.Model

    """
    # We're going to import the Matlab model and convert it to COBApy
    mat_data = loadmat(path)

    if variable_name:
        mat_model = mat_data[variable_name][0, 0]
    else:
        meta_vars = {"__globals__", "__header__", "__version__"}
        possible_keys = sorted(i for i in mat_data if i not in meta_vars)

        if len(possible_keys) == 1:
            variable_name = possible_keys[0]
        else:
            raise Exception('Need to specify a variable name if several ' + \
                            'variables are contained int the mat file')

        mat_model = mat_data[variable_name][0,0]


    cobra_model = Model(mat_model['description'][0])

    cobra_model.description = mat_model['description'][0]
    ## METABOLITES
    # In the Matlab model, the corresponding components are :
    # * mets = Identifiers
    # * metNames = Names
    # * metFormulas = formulas
    # * metCompSymbol = compartments

    metabolites = [
        Metabolite(
            mat_model['mets'][i,0][0],
            formula=mat_model['metFormulas'][i,0][0],
            name=mat_model['metNames'][i,0][0],
            compartment=mat_model['metCompSymbol'][i,0][0]
        ) for i in range(len(mat_model['metNames']))
    ]

    # Get the metSEEDID
    seed_id = mat_model['metSEEDID']
    for i,_ in enumerate(metabolites):
        metabolites[i].annotation = {"seed_id":seed_id[i,0][0]}

    ## REACTIONS
    # In the Matlab model, the corresponding components are :
    # * rxns = Names
    # * rev = Reversibility (not used, see https://cobrapy.readthedocs.io/en/0.5.11/faq.html#How-do-I-change-the-reversibility-of-a-Reaction?)
    # * lb = Lower bounds
    # * ub = Upper bounds
    # * subSystems = subsystem names
    # * S : Reactions matrix
    #       1 line = 1 metabolite
    #       1 column = 1 reaction
    # * c : Objective coefficient
    # * genes : Genes names
    # * rules : Genes rules

    # Some utilities for gene generation rules (convert the IDs to the names of the genes)
    gene_pattern = re.compile(r'x\([0-9]+\)')

    def gene_id_to_name(match):
        id_ = int(match.group()[2:-1])
        return mat_model['genes'][id_,0][0]

    # Add each reaction
    for i in range(mat_model['S'].shape[1]):
        # Name of the reaction
        reaction = Reaction(str(mat_model['rxns'][i,0][0]))

        # Add the reaction to the model
        cobra_model.add_reaction(reaction)

        # NOTE : The str() conversion above is needed, otherwise the CPLEX solver
        # does not work : "Invalid matrix input type --"

        # Reaction description
        # reaction.name not set

        # Subsystem (only if set in the Matlab model)
        if(len(mat_model['subSystems'][i,0])):
            reaction.subsystem = mat_model['subSystems'][i,0][0]

        # Lower bound
        reaction.lower_bound = float(mat_model['lb'][i,0])
        # Upper bound
        reaction.upper_bound = float(mat_model['ub'][i,0])

        # Objective coefficient
        reaction.objective_coefficient = float(mat_model['c'][i,0])

        # Metabolites
        react_mets = {}
        # Iterate over each metabolite and see if it is part of the reaction
        # (stoechiomectric coefficient not equal to 0)
        for j,_ in enumerate(metabolites):
            if(mat_model['S'][j,i] != 0):
                react_mets[metabolites[j]] = mat_model['S'][j,i]

        reaction.add_metabolites(react_mets)

        # Genes
        try:
            if len(mat_model['rules'][i,0]):
                rule = mat_model['rules'][i,0][0]
                # Call the regex magic to convert IDs to gene names
                rule = gene_pattern.sub(gene_id_to_name, rule)

                # Add the data to the reaction
                reaction.gene_reaction_rule = rule
        except ValueError:
            pass

    Compartments = dict()
    CompartmentDB = mat_model['CompartmentData'][0]

    num_comp = len(CompartmentDB['pH'][0][0])
    comps = [{} for i in range(num_comp)]

    for i in range(num_comp):
        comp = comps[i]
        comp['pH'] = CompartmentDB['pH'][0][0][i]
        comp['ionicStr'] = CompartmentDB['ionicStr'][0][0][i]
        comp['symbol'] = CompartmentDB['compSymbolList'][0][0,i][0]
        comp['name'] = CompartmentDB['compNameList'][0][0,i][0]
        comp['c_max'] = CompartmentDB['compMaxConc'][0][0][i]
        comp['c_min'] = CompartmentDB['compMinConc'][0][0][i]

        # Register our compartment by its name
        if comp['symbol'] in Compartments:
            raise Exception("Duplicate compartment ID : " + comp['symbol'])

        Compartments[comp['symbol']] = comp

    # We need to iterate first once to set the names of the compartments, so that we have the keys for our dictionnary...
    for i in range(num_comp):
        comp = comps[i]
        comp['membranePot'] = {}
        for j in range(num_comp):
            comp['membranePot'][comps[j]['symbol']] = CompartmentDB['membranePot'][0][i,j]

    cobra_model.compartments = Compartments

    return cobra_model


def load_thermoDB(path):
    """ Load a thermodynamic database

    :param string path: The path of the file to load
    :returns: The thermodynamic database
    :rtype: dict

    """
    with open(path, 'rb') as file:
        ReactionDB = pickle.loads(zlib.decompress(file.read()))

    return ReactionDB


def printLP(model):
    """ Print the LP file corresponding to the model

    :param cobra.core.model.Model model: The model to output the LP file for

    :returns: The content of the LP file
    :rtype: str

    Usually, you pass the result of this function to :py:func:`file.write` to
    write it on disk. If you prefer, you can use :func:`pytfa.io.writeLP` to
    write the result directly to a file.

    """
    # It's much more efficient to allocate the file in memory then write it as
    # a whole to the disk, rather than making multiple disk calls
    res = ''

    # Write the problem name
    res += '\\Problem name: {}_LP\n\nMaximize\n obj: '.format(model.description)

    # Write the objective
    for rxn in model.reactions:
        if rxn.objective_coefficient != 0:
            if rxn.objective_coefficient != 1:
                res += rxn.objective_coefficient + ' '
            res += rxn.id

    # Constraints
    res += '\nSubject To\n'

    for cons in model.constraints:
        # Name of the constraint :
        res += cons.name + ': '

        # Write the lower bound if applicable...
        if cons.lb != cons.ub:
            if cons.lb != None:
                res += str(cons.lb)  + ' < '

        # Write the bound
        res += str(cons.expression)

        # Write the upper bound
        if cons.lb == cons.ub:
            res += ' = ' + str(cons.ub)
        elif cons.ub != None:
            res += ' < ' + str(cons.ub)

        # Next line
        res += '\n'

    # Variables
    res += 'Bounds\n'

    for var in model.variables:
        # optlang already does the hard job for us, yay !
        res += str(var) + '\n'

    # Binrary constraints
    res += 'Binaries\n'

    # Number of variables on the current line we're writing
    count = 1

    for var in model.variables:
        if var.type=='binary':
            res += var.name + '\t'
            count += 1
            # Print at most 7 variables per line
            if count == 7:
                res += '\n'
                count = 1

    # Done !
    res += '\nEnd'

    return res


def writeLP(model, path = None):
    """ Write the LP file of the specified model to the file indicated by path.

    :param cobra.core.model.Model model: The COBRApy model to write the LP file
        for
    :param string path: `Optional` The path of the file to be written. If not
        specified, the name of the COBRApy model will be used.


    """
    if not path:
        path = model.description + '.lp'

    with open(path, 'w') as file:
        file.write(printLP(model))


