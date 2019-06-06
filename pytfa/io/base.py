"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Input/Output tools to import or export pytfa models


"""

import pickle
import zlib
import numpy as np
import re

from cobra import Model, Reaction, Metabolite
from cobra.io import load_matlab_model
from scipy.io import loadmat, savemat

from ..utils.numerics import BIGM_DG

from warnings import warn

try:
    from scipy.sparse import dok_matrix, lil_matrix
except ImportError:
    dok_matrix, lil_matrix = None, None


def import_matlab_model(path, variable_name=None):
    """ Convert at matlab cobra_model to a pyTFA cobra_model, with Thermodynamic values

    :param variable_name:
    :param string path: The path of the file to import

    :returns: The converted cobra_model
    :rtype: cobra.thermo.model.Model

    """
    # We're going to import the Matlab cobra_model and convert it to COBApy
    mat_data = loadmat(path)

    if variable_name:
        mat_model = mat_data[variable_name][0, 0]
    else:
        meta_vars = {"__globals__", "__header__", "__version__"}
        possible_keys = sorted(i for i in mat_data if i not in meta_vars)

        if len(possible_keys) == 1:
            variable_name = possible_keys[0]
        else:
            raise Exception(
                'Need to specify a variable name if several ' + 'variables are contained int the mat file')

        mat_model = mat_data[variable_name][0, 0]

    # cobra_model = Model(mat_model['description'][0])
    cobra_model = load_matlab_model(path)

    cobra_model.description = mat_model['description'][0]
    cobra_model.name = mat_model['description'][0]

    metabolites = cobra_model.metabolites

    ## METABOLITES
    # In the Matlab cobra_model, the corresponding components are :
    # * mets = Identifiers
    # * metNames = Names
    # * metFormulas = formulas
    # * metCompSymbol = compartments

    # def read_mat_model(mat_struct, field_name, index):
    #     try:
    #         return mat_struct[field_name][index,0][0]
    #     except IndexError:
    #         return None

    # metabolites = [Metabolite( read_mat_model(mat_model,'mets',i),
    #     formula=read_mat_model(mat_model,'metFormulas',i),
    #     name=read_mat_model(mat_model,'metNames',i),
    #     compartment=read_mat_model(mat_model,'metCompSymbol',i))
    #                for i in range(len(mat_model['metNames']))]

    # Get the metSEEDID
    seed_id = mat_model['metSEEDID']
    for i, _ in enumerate(metabolites):
        metabolites[i].annotation = {"seed_id": seed_id[i, 0][0]}

    # ## REACTIONS
    # # In the Matlab cobra_model, the corresponding components are :
    # # * rxns = Names
    # # * rev = Reversibility (not used, see https://cobrapy.readthedocs.io/en/0.5.11/faq.html#How-do-I-change-the-reversibility-of-a-Reaction?)
    # # * lb = Lower bounds
    # # * ub = Upper bounds
    # # * subSystems = subsystem names
    # # * S : Reactions matrix
    # #       1 line = 1 metabolite
    # #       1 column = 1 reaction
    # # * c : Objective coefficient
    # # * genes : Genes names
    # # * rules : Genes rules

    # # Some utilities for gene generation rules (convert the IDs to the names of the genes)
    # gene_pattern = re.compile(r'x\([0-9]+\)')

    # def gene_id_to_name(match):
    #     id_ = int(match.group()[2:-1])
    #     # /!\ These are indexed from 1, while python indexes from 0
    #     return mat_model['genes'][id_ - 1, 0][0]

    # # Add each reaction
    # for i in range(mat_model['S'].shape[1]):
    #     # Name of the reaction
    #     reaction = Reaction(str(mat_model['rxns'][i, 0][0]))

    #     # Add the reaction to the cobra_model
    #     cobra_model.add_reactions([reaction])

    #     # NOTE : The str() conversion above is needed, otherwise the CPLEX solver
    #     # does not work : "Invalid matrix input type --"

    #     # Reaction description
    #     # reaction.name not set

    #     # Subsystem (only if set in the Matlab cobra_model)
    #     if (len(mat_model['subSystems'][i, 0])):
    #         reaction.subsystem = mat_model['subSystems'][i, 0][0]

    #     # Lower bound
    #     reaction.lower_bound = float(mat_model['lb'][i, 0])
    #     # Upper bound
    #     reaction.upper_bound = float(mat_model['ub'][i, 0])

    #     # Objective coefficient
    #     reaction.objective_coefficient = float(mat_model['c'][i, 0])

    #     # Metabolites
    #     react_mets = {}
    #     # Iterate over each metabolite and see if it is part of the reaction
    #     # (stoechiomectric coefficient not equal to 0)
    #     for j, _ in enumerate(metabolites):
    #         if (mat_model['S'][j, i] != 0):
    #             react_mets[metabolites[j]] = mat_model['S'][j, i]

    #     reaction.add_metabolites(react_mets)

    #     # Genes
    #     try:
    #         if len(mat_model['rules'][i, 0]):
    #             rule = mat_model['rules'][i, 0][0]
    #             # Call the regex magic to convert IDs to gene names
    #             rule = gene_pattern.sub(gene_id_to_name, rule)

    #             # Add the data to the reaction
    #             reaction.gene_reaction_rule = rule
    #     except ValueError:
    #         pass
    #     except IndexError:
    #         # The gene number is higher than the length of the gene list
    #         warn('The gene reaction rule {} appears to be misaligned with '
    #              'the gene list'.format(rule))


    Compartments = dict()
    CompartmentDB = mat_model['CompartmentData'][0]

    num_comp = len(CompartmentDB['pH'][0][0])
    comps = [{} for i in range(num_comp)]

    for i in range(num_comp):
        comp = comps[i]
        comp['pH'] = CompartmentDB['pH'][0][0][i]
        comp['ionicStr'] = CompartmentDB['ionicStr'][0][0][i]
        comp['symbol'] = CompartmentDB['compSymbolList'][0][0, i][0]
        comp['name'] = CompartmentDB['compNameList'][0][0, i][0]
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
            comp['membranePot'][comps[j]['symbol']] = \
            CompartmentDB['membranePot'][0][i, j]

    cobra_model.compartments = Compartments

    return cobra_model


def write_matlab_model(tmodel, path, varname='tmodel'):
    """
    Writes the Thermo Model to a Matlab-compatible structure

    :param varname:
    :param tmodel:
    :param path:
    :return: None
    """

    from cobra.io.mat import create_mat_dict

    mat = create_mat_dict(tmodel)
    mat.update(create_thermo_dict(tmodel))
    mat.update(create_problem_dict(tmodel))
    savemat(path, {varname: mat}, appendmat=True,
                     oned_as="column")


def create_thermo_dict(tmodel):
    """
    Dumps the thermodynamic information in a mat-compatible dictionary
    (similar to the output of cobra.io.mat.create_mat_dict)

    :param tmodel: pytfa.thermo.tmodel.ThermoModel
    :return: dict object
    """

    mat = {}

    met_map = {
        'metSEEDID':('id',''),
        'metDeltaGFstd':('deltaGf_std',BIGM_DG),
        'metDeltaGFerr':('deltaGf_err',BIGM_DG),
        'metMass':('mass',BIGM_DG),
        'metCharge':('charge_std',BIGM_DG),
        'metDelGFtr':('deltaGf_tr',BIGM_DG),
        'metCompSymbol':('compartment','')
    }

    for column,(key,default_value) in met_map.items():
        the_data = np.full(len(tmodel.metabolites),default_value)

        for e,the_met in enumerate(tmodel.metabolites):
            try:
                the_value = the_met.thermo[key]
                the_data[e] = the_value
            except KeyError:
                continue

        if column == 'metSEEDID':
            mat[column] = np.array(the_data, dtype=np.object)
        else:
            mat[column] = np.array(the_data)


    rxn_map = {
        'rxnDeltaGR':('DeltaGR',BIGM_DG),
        'rxnDeltaGRerr':('DeltaGRerr',BIGM_DG),
        'rxnThermo':('computed',None),
        'isTrans':('isTrans',None),
    }

    for column,(key,default_value) in rxn_map.items():
        the_data = np.full(len(tmodel.reactions),default_value)

        for e,the_rxn in enumerate(tmodel.reactions):
            try:
                the_value = the_rxn.thermo[key]
                the_data[e] = the_value
            except KeyError:
                continue

        mat[column] = np.array(the_data)

    # Adding compartment data
    CompartmentDB = {}

    compartments = [v for v in tmodel.compartments.values()]
    CompartmentDB['pH'] = np.array([x['pH'] for x in compartments])
    CompartmentDB['ionicStr'] = np.array([x['ionicStr'] for x in compartments])
    CompartmentDB['compMaxConc'] = np.array([x['c_max'] for x in compartments])
    CompartmentDB['compMinConc'] = np.array([x['c_min'] for x in compartments])

    #Write symbols and names in collumn cell arrays
    CompartmentDB['compSymbolList'] = np.zeros((1, len(compartments)), dtype=np.object)
    CompartmentDB['compNameList'] =  np.zeros((1, len(compartments)), dtype=np.object)

    mat_to_python_string = [('compSymbolList', 'symbol'),
                            ('compNameList', 'name')]

    for i, _ in enumerate(compartments):
        for mat_name, py_name in mat_to_python_string:
            CompartmentDB[mat_name][0, i] = compartments[i][py_name]


    # The membrane potential is an NxN matrix in the matlab format
    membrane_pot = np.zeros((len(compartments),len(compartments)))
    for i, _ in enumerate(compartments):
        for j, _ in enumerate(compartments):
            symbol_j = compartments[j]['symbol']
            membrane_pot[i,j] = compartments[i]['membranePot'][symbol_j]

    CompartmentDB['membranePot'] = membrane_pot

    mat['CompartmentData'] = CompartmentDB

    return mat


def varnames2matlab(name, tmodel):
    """
    Transforms reaction variable pairs from `('ACALD','ACALD_reverse_xxxxx')` to
    `('F_ACALD','B_ACALD')`   if it is a reaction, else leaves is as is

    :return:
    """

    reverse_regex = re.compile(r'(.+_reverse)_[a-f0-9]{5}')

    new_name = name

    if new_name in tmodel.reactions:
        new_name = 'F_' + new_name
    else:
        test = reverse_regex.match(new_name)
        if test:
            new_name = 'R_' + test.groups()[0]

    return new_name



def create_problem_dict(tmodel):
    """
    Dumps the the MILP formulation for TFA in a mat-compatible dictionary
    (similar to the output of cobra.io.mat.create_mat_dict)

    :param tmodel: pytfa.thermo.tmodel.ThermoModel
    :ret
    """

    vartype_map = {
        'continuous':'C',
        'binary':'B',
        'integer':'B' # FIXME Waiting for optlang fix on binary vars
    }

    mat= {}

    mat['A'] = create_generalized_matrix(tmodel)

    # Variables

    mat['var_lb'] = np.array([x.lb for x in tmodel.variables]) * 1.
    mat['var_ub'] = np.array([x.ub for x in tmodel.variables]) * 1.
    vname = np.full(len(tmodel.variables),'', dtype=np.object)
    vtype = np.full(len(tmodel.variables),'', dtype=np.object)

    for e, this_var in enumerate(tmodel.variables):
        this_name = this_var.name
        new_name = varnames2matlab(this_name, tmodel)

        vname[e] = new_name
        vtype[e] = vartype_map[this_var.type]

    mat['varNames'] = vname
    mat['vartypes'] = vtype

    obj = tmodel.objective
    mat['objtype'] = -1 if tmodel.objective.direction.startswith('max') else 1
    mat['f'] = np.full(len(tmodel.variables), 0, dtype = np.double) * 1.0

    for this_var,this_coeff in obj.get_linear_coefficients(obj.variables).items():
        matlab_name = varnames2matlab(this_var.name, tmodel)
        mat['f'][np.where(mat['varNames']==matlab_name)[0][0]] = this_coeff

    # Constraints

    rhs = np.empty(len(tmodel.constraints))
    ctype = np.full(len(tmodel.constraints),'', dtype=np.object)
    cname = np.full(len(tmodel.constraints),'', dtype=np.object)

    for e,this_cons in enumerate(tmodel.constraints):
        if   this_cons.lb is None and this_cons.ub is not None:
            this_type = '<'
            this_rhs = this_cons.ub
        elif this_cons.ub is None and this_cons.lb is not None:
            this_type = '>'
            this_rhs = this_cons.lb
        elif this_cons.ub is not None and this_cons.lb is not None:
            if this_cons.lb == this_cons.ub:
                this_type = '='
                this_rhs = this_cons.ub
            else:
                raise(Exception('Matlab formulation does not support '
                                'constraints under the form: lb<x<ub'))
        else:
            raise(ValueError('Constraint type not recognized'))

        rhs[e] = this_rhs
        ctype[e] = this_type
        cname[e] = this_var.name

    mat['rhs'] = rhs *1
    mat['constraintType'] = ctype
    mat['constraintNames'] = cname

    return mat


def create_generalized_matrix(tmodel, array_type = 'dense'):
    """
    Returns the generalized stoichiomatric matrix used for TFA

    :param array_type:
    :param tmodel: pytfa.ThermoModel

    :returns: matrix.
    """

    if array_type not in ('DataFrame', 'dense') and not dok_matrix:
        raise ValueError('Sparse matrices require scipy')

    dtype = np.float64

    array_constructor = {'dense': np.zeros, 'dok': dok_matrix,
        'lil': lil_matrix, 'DataFrame': np.zeros, }

    n_constraints = len(tmodel.constraints)
    n_variables = len(tmodel.variables)
    array = array_constructor[array_type]((n_constraints, n_variables),
                                          dtype=dtype)

    c_ind = {x:e for e,x in enumerate(tmodel.constraints)}
    v_ind = {x:e for e,x in enumerate(tmodel.variables)}

    for this_cons in tmodel.constraints:
        var_coeff_dict = this_cons.get_linear_coefficients(this_cons.variables)

        for this_var,coeff in var_coeff_dict.items():
            array[c_ind[this_cons], v_ind[this_var]] = coeff

    if array_type == 'DataFrame':
        metabolite_ids = [met.id for met in tmodel.constraints]
        reaction_ids = [rxn.id for rxn in tmodel.variables]
        return pd.DataFrame(array, index=metabolite_ids, columns=reaction_ids)

    else:
        return array


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
    """ Print the LP file corresponding to the cobra_model

    :param cobra.thermo.model.Model model: The cobra_model to output the LP file for

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
            if cons.lb is not None:
                res += str(cons.lb) + ' < '

        # Write the bound
        res += str(cons.expression)

        # Write the upper bound
        if cons.lb == cons.ub:
            res += ' = ' + str(cons.ub)
        elif cons.ub is not None:
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
        # FIXME using the integer trick to be able to constraint binary variables
        # if var.type=='binary':
        if var.type in ['binary', 'integer']:
            res += var.name + '\t'
            count += 1
            # Print at most 7 variables per line
            if count == 7:
                res += '\n'
                count = 1

    # Done !
    res += '\nEnd'

    return res


def writeLP(model, path=None):
    """ Write the LP file of the specified cobra_model to the file indicated by path.

    :param cobra.thermo.model.Model model: The COBRApy cobra_model to write the LP file
        for
    :param string path: `Optional` The path of the file to be written. If not
        specified, the name of the COBRApy cobra_model will be used.


    """
    if not path:
        path = model.description + '.lp'

    with open(path, 'w') as file:
        file.write(printLP(model))
