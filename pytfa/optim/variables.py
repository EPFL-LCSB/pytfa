# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Variable declarations

"""
from distutils.version import LooseVersion
from optlang import __version__ as OPTLANG_VER

from ..utils.str import camel2underscores

from warnings import warn


op_replace_dict = {
    ' + ': '_ADD_',
    ' - ': '_SUB_',
    ' * ': '_MUL_',
    ' / ': '_DIV_',
}

###################################################
###                  VARIABLES                  ###
###################################################


class GenericVariable:
    """
    Class to represent a generic variable. The purpose is that the interface
    is instantiated on initialization, to follow the type of interface used
    by the problem, and avoid incompatibilities in optlang

    Attributes:

        :id: Used for DictList comprehension. Usually points back at a
        enzyme or reaction id for ease of linking. Should be unique given
        a variable type.
        :name: Should be a concatenation of the id and a prefix that is
        specific to the variable type. will be used to address the variable at
        the solver level, and hence should be unique in the whole cobra_model
        :cobra_model: the cobra_model hook.
        :variable: links directly to the cobra_model representation of tbe variable
    """
    prefix = ''

    @property
    def __attrname__(self):
        """
        Name the attribute the instances will have
        Example: GenericVariable -> generic_variable
        :return:
        """
        return camel2underscores(self.__class__.__name__)

    def __init__(self, id_='', model=None, hook=None, queue=False, scaling_factor=1,
                 **kwargs):
        """

        :param id_: will be used to identify the variable
            (name will be a concat of this and a prefix)
        :param model: the cobra.Model object
        :param queue: whether or not to queue the variable for update object
        :param scaling_factor: scaling factor used in self.scaled, useful for
                                adimensionalisation of constraints
        :param kwargs: stuff you want to pass to the variable constructor
        """
        self.hook = hook
        self._id = id_
        self._model = model
        self.kwargs = kwargs
        self._name = self.make_name()
        self.get_interface(queue)
        self.prefix = ''
        self._scaling_factor = scaling_factor

    def get_interface(self, queue):
        """
        Called upon completion of __init__, initializes the value of self.var,
        which is returned upon call, and stores the actual interfaced variable.

        :return: instance of Variable from the problem
        """

        if not self.name in self.model.variables:
            variable = self.model.problem.Variable(name = self.name, **self.kwargs)
            if not queue:
                self.model.add_cons_vars(variable, sloppy=self.model.sloppy)
            else:
                self.model._var_queue.append(variable)
        else:
            self.variable = self.model.variables.get(self.name)

    def make_name(self):
        """
        Needs to be overridden by the subclass, concats the id with a
         prefix

        :return: None
        """
        return self.prefix + self.id

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def id(self):
        """
        for cobra.thermo.DictList compatibility
        :return:
        """
        return self._id

    @property
    def variable(self):
        return self.model.variables[self.name]

    @variable.setter
    def variable(self,value):
        self.model.variables[self.name] = value

    @property
    def scaling_factor(self):
        return self._scaling_factor

    @property
    def unscaled(self):
        """
        If the scaling factor of quantity X is a, it is represented by
        the variable X_hat = X/a.
        This returns X = a.X_hat
        Useful for nondimensionalisation of variables and constraints

        :return: The variable divided by its scaling factor
        """
        return self.scaling_factor * self

    @property
    def value(self):
        try:
            return self.model.solution.x_dict[self.name]
        except AttributeError:
            warn('Model need to be optimized to get a value for this variable')

    @property
    def unscaled_value(self):
        try:
            return self.scaling_factor * self.value
        except AttributeError:
            warn('Model need to be optimized to get a value for this variable')

    @property
    def model(self):
        return self._model
    #
    # @variable.setter
    # def variable(self, value):
    #     self._variable = value

    @property
    def type(self):
        return self.variable.type

    def test_consistency(self, other):
        """
        Tests whether a candidate to an operation is of the right type and is
        from the same problem

        :param other: an object
        :return: None
        """

        assert (isinstance(other, GenericVariable))
        # assert (other.problem == self.problem)

    def get_operand(self, other):
        """
        For operations, choose if the operand is a GenericVariable, in which
        we return its optlang variable, or something else (presumably a numeric)
        and we let optlang decide what to do

        :param other:
        :return:
        """
        if isinstance(other, GenericVariable):
            return other.variable
        else:  # Let optlang decide what to do
            return other

    #########################################################################
    # We redefine all the following operations to return the equivalent     #
    # expression one would get by applying the said operation to the        #
    # self.variable attribute                                               #
    #########################################################################

    def __add__(self, other):
        """
        Adding either two variables together or a variable and a numeric
        results in a new variable
        :param other:
        :return: a new Generic Variable
        """
        operand = self.get_operand(other)

        new_variable = self.variable + operand

        return self.make_result(new_variable)

    def __radd__(self, other):
        """
        Take priority on symmetric arithmetic operation
        :param other:
        :return:
        """

        return self.__add__(other)

    def __sub__(self, other):
        """
        Substracting either two variables together or a variable and a numeric
        results in a new variable
        :param other:
        :return: a new Generic Variable
        """
        operand = self.get_operand(other)

        new_variable = self.variable - operand

        return self.make_result(new_variable)

    def __rsub__(self, other):
        """
        Take priority on symmetric arithmetic operation
        :param other:
        :return:
        """

        operand = self.get_operand(other)

        new_variable = operand - self.variable

        return self.make_result(new_variable)

    def __mul__(self, other):
        """
        Multiplying either two variables together or a variable and a numeric
        results in a new variable
        :param other:
        :return: a new Generic Variable
        """
        operand = self.get_operand(other)

        new_variable = self.variable * operand

        return self.make_result(new_variable)

    def __rmul__(self, other):
        """
        Take priority on symmetric arithmetic operation
        :param other:
        :return:
        """
        return self.__mul__(other)

    def __truediv__(self, other):
        """
        Dividing either two variables together or a variable and a numeric
        results in a new variable
        :param other:
        :return: a new Generic Variable
        """
        operand = self.get_operand(other)

        new_variable = self.variable / operand

        return self.make_result(new_variable)

    def __rtruediv__(self, other):
        """
        Take priority on symmetric arithmetic operation
        :param other:
        :return:
        """

        operand = self.get_operand(other)

        new_variable = operand / self.variable

        return self.make_result(new_variable)


    def make_result(self, new_variable):
        """
        Returns a Sympy expression
        :param new_variable:
        :return:
        """

        # if isinstance(new_variable, self.cobra_model.problem.Variable):
        #     new_name = new_variable.name
        # else: #It is an expression
        #     new_name  = new_variable.__str__()
        #     for k,v in op_replace_dict.items():
        #         new_name = new_name.replace(k,v)

        # result = GenericVariable(new_name, self.cobra_model)
        # result.variable = new_variable

        # return result
        return new_variable

    def __repr__(self):
        return self.variable.__repr__()

def get_binary_type():
    """
    FIX : We enforce type to be integer instead of binary, else optlang does
    not allow to set the binary variable bounds to anything other than (0,1)
    You might want to set it at (0,0) to enforce directionality for example
    """
    if LooseVersion(OPTLANG_VER) < LooseVersion('1.2.3'):
        return 'integer'
    else:
        return 'binary'

class ModelVariable(GenericVariable):
    """
    Class to represent a variable attached to the model
    """

    def __init__(self, model, id_, **kwargs):
        if not 'lb' in kwargs:
            kwargs['lb'] = 0
        GenericVariable.__init__(self,
                                 id_= id_,
                                 model=model,
                                 hook=model,
                                 **kwargs)

class BinaryVariable(GenericVariable):
    """
    Class to represent a generic binary variable
    """

    def __init__(self, id_, model, **kwargs):

        if not 'lb' in kwargs:
            kwargs['lb'] = 0
        if not 'ub' in kwargs:
            kwargs['ub'] = 1

        GenericVariable.__init__(self,
                                 id_,
                                 model = model,
                                 type=get_binary_type(),
                                 **kwargs)

    prefix = 'B_'

class ReactionVariable(GenericVariable):
    """
    Class to represent a variable attached to a reaction
    """

    def __init__(self, reaction, **kwargs):
        model = reaction.model

        GenericVariable.__init__(self,
                                 model=model,
                                 hook=reaction,
                                 **kwargs)

    @property
    def reaction(self):
        return self.hook

    @property
    def id(self):
        return self.reaction.id

    @property
    def model(self):
        return self.reaction.model

class MetaboliteVariable(GenericVariable):
    """
    Class to represent a variable attached to a enzyme
    """

    def __init__(self, metabolite, **kwargs):
        model = metabolite.model

        GenericVariable.__init__(self,
                                 model=model,
                                 hook=metabolite,
                                 **kwargs)

    @property
    def metabolite(self):
        return self.hook

    @property
    def id(self):
        return self.metabolite.id

    @property
    def model(self):
        return self.metabolite.model


class ForwardUseVariable(ReactionVariable, BinaryVariable):
    """
    Class to represent a forward use variable, a type of binary variable used to
    enforce forward directionality in reaction net fluxes
    """

    def __init__(self, reaction, **kwargs):
        ReactionVariable.__init__(self, reaction,
                                  type=get_binary_type(),
                                  **kwargs)

    prefix = 'FU_'

class BackwardUseVariable(ReactionVariable, BinaryVariable):
    """
    Class to represent a backward use variable, a type of binary variable used
    to enforce backward directionality in reaction net fluxes
    """

    def __init__(self, reaction, **kwargs):
        if not 'lb' in kwargs:
            kwargs['lb'] = 0
        if not 'ub' in kwargs:
            kwargs['ub'] = 1

        ReactionVariable.__init__(self, reaction,
                                  type=get_binary_type(),
                                  **kwargs)

    prefix = 'BU_'


class LogConcentration(MetaboliteVariable):
    """
    Class to represent a log concentration of a enzyme
    """

    prefix = 'LC_'


class DeltaGErr(ReactionVariable):
    """
    Class to represent a DeltaGErr
    """

    prefix = 'DGE_'


class DeltaG(ReactionVariable):
    """
    Class to represent a DeltaG
    """

    prefix = 'DG_'


class DeltaGstd(ReactionVariable):
    """
    Class to represent a DeltaG^o (naught) - standard conditions
    """

    prefix = 'DGo_'

class ThermoDisplacement(ReactionVariable):
    """
    Class to represent the thermodynamic displacement of a reaction
    \Gamma = -\DeltaG/RT
    """

    prefix = 'LnGamma_'

class PosSlackVariable(ReactionVariable):
    """
    Class to represent a positive slack variable for relaxation problems
    """

    def __init__(self,reaction,**kwargs):
        ReactionVariable.__init__(self, reaction, **kwargs)

    prefix = 'PosSlack_'


class NegSlackVariable(ReactionVariable):
    """
    Class to represent a negative slack variable for relaxation problems
    """

    def __init__(self,reaction,**kwargs):
        ReactionVariable.__init__(self, reaction, **kwargs)

    prefix = 'NegSlack_'

class PosSlackLC(MetaboliteVariable):

    prefix = 'PosSlackLC_'

class NegSlackLC(MetaboliteVariable):

    prefix = 'NegSlackLC_'

