"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Sampling wrappers for pytfa models


"""

import numpy as np
from sympy.core.singleton import S
from time import  time
from  cobra.sampling import OptGPSampler, ACHRSampler, HRSampler,\
                                            shared_np_array
from optlang.interface import OPTIMAL

class GeneralizedHRSampler(HRSampler):

    def __init__(self, model, thinning,  nproj=None, seed=None):
        """
        Adapted from cobra.flux_analysis.sampling.py
        _________________________________________

        Initialize a new sampler object.
        """
        # This currently has to be done to reset the solver basis which is
        # required to get deterministic warmup point generation
        # (in turn required for a working `seed` arg)
        HRSampler.__init__(self, model, thinning, seed=seed)
        
        if model.solver.is_integer:
            raise TypeError("sampling does not work with integer problems :(")
        self.model = model.copy()

        self.thinning = thinning

        if nproj is None:
            self.nproj = int(min(len(self.model.variables)**3, 1e6))
        else:
            self.nproj = nproj
        self.n_samples = 0
        self.retries = 0
        
        # Careful, double underscore names are mangled to the class name
        self.problem = self._HRSampler__build_problem()

        # Set up a map from reaction -> forward/reverse variable
        var_idx = {v: idx for idx, v in enumerate(self.model.variables)}
        self.var_idx = np.array(
            [var_idx[v] for v in self.model.variables])
        self.warmup = None
        if seed is None:
            self._seed = int(time())
        else:
            self._seed = seed
        # Avoid overflow
        self._seed = self._seed % np.iinfo(np.int32).max

    def generate_fva_warmup(self):
        """
        Adapted from cobra.flux_analysis.sampling.py
        __________________________________________

        Generate the warmup points for the sampler.

        Generates warmup points by setting each flux as the sole objective
        and minimizing/maximizing it. Also caches the projection of the
        warmup points into the nullspace for non-homogeneous problems (only
        if necessary).
        """
        self.n_warmup = 0
        idx = np.hstack([self.var_idx])
        self.warmup = np.zeros((len(idx), len(self.model.variables)))
        self.model.objective = S.Zero
        self.model.objective.direction = "max"
        variables = self.model.variables
        for i in idx:
            # Omit fixed reactions
            if self.problem.variable_fixed[i]:
                self.model.logger.info("skipping fixed variable %s" % variables[i].name)
                continue
            self.model.objective.set_linear_coefficients({variables[i]: 1})
            self.model.slim_optimize()
            if not self.model.solver.status == OPTIMAL:
                self.model.logger.info(
                    "can not maximize variable %s, skipping it" % variables[
                        i].name)
                continue
            primals = self.model.solver.primal_values
            sol = [primals[v.name] for v in self.model.variables]
            self.warmup[self.n_warmup,] = sol
            self.n_warmup += 1
            # revert objective
            self.model.objective.set_linear_coefficients({variables[i]: 0})
        # Shrink warmup points to measure
        self.warmup = shared_np_array((self.n_warmup, len(variables)),
                                      self.warmup[0:self.n_warmup, ])


# Next, we redefine the analysis class as both inheriting from the
# GeneralizedHRSampler, and then the original classes. This will overwrite the
# inherited methods from HRSampler, while keeping the lower level ones from the
# samplers

class GeneralizedACHRSampler(GeneralizedHRSampler,ACHRSampler):
    def __init__(self, model, thinning=100, seed=None):
        """
        Adapted from cobra.flux_analysis.analysis
        __________________________________________
        Initialize a new ACHRSampler."""
        GeneralizedHRSampler.__init__(self, model, thinning, seed=seed)
        self.generate_fva_warmup()
        self.prev = self.center = self.warmup.mean(axis=0)
        np.random.seed(self._seed)

class GeneralizedOptGPSampler(GeneralizedHRSampler, OptGPSampler):
    def __init__(self, model, processes, thinning=100, seed=None):
        """
        Adapted from cobra.flux_analysis.sampling.py
        __________________________________________
        Initialize a new OptGPSampler."""
        GeneralizedHRSampler.__init__(self, model, thinning, seed=seed)
        self.generate_fva_warmup()
        self.processes = processes

        # This maps our saved center into shared memory,
        # meaning they are synchronized across processes
        self.center = shared_np_array((len(self.model.variables), ),
                                      self.warmup.mean(axis=0))


def sample(model, n, method="optgp", thinning=100, processes=1, seed=None):
    """
    Sample valid flux distributions from a thermo cobra_model.

    Function adapted from cobra.flux_analysis.sample to display all solver
    variables

    **Documentation adapted from cobra.flux_analysis.sample**

    1. 'optgp' (default) which uses the OptGPSampler that supports parallel
        analysis [1]_. Requires large numbers of samples to be performant
        (n < 1000). For smaller samples 'achr' might be better suited.

    or

    2. 'achr' which uses artificial centering hit-and-run. This is a single
       process method with good convergence [2]_.

    Parameters
    ----------
    model : pytfa.core.ThermoModel
        The cobra_model from which to sample variables.
    n : int
        The number of samples to obtain. When using 'optgp' this must be a
        multiple of `processes`, otherwise a larger number of samples will be
        returned.
    method : str, optional
        The analysis algorithm to use.
    thinning : int, optional
        The thinning factor of the generated analysis chain. A thinning of 10
        means samples are returned every 10 steps. Defaults to 100 which in
        benchmarks gives approximately uncorrelated samples. If set to one
        will return all iterates.
    processes : int, optional
        Only used for 'optgp'. The number of processes used to generate
        samples.
    seed : positive integer, optional
        The random number seed to be used. Initialized to current time stamp
        if None.

    Returns
    -------
    pandas.DataFrame
        The generated flux samples. Each row corresponds to a sample of the
        fluxes and the columns are the reactions.

    Notes
    -----
    The samplers have a correction method to ensure equality feasibility for
    long-running chains, however this will only work for homogeneous models,
    meaning models with no non-zero fixed variables or constraints (
    right-hand side of the equalities are zero).

    References
    ----------
    .. [1] Megchelenbrink W, Huynen M, Marchiori E (2014)
       optGpSampler: An Improved Tool for Uniformly Sampling the Solution-Space
       of Genome-Scale Metabolic Networks.
       PLoS ONE 9(2): e86587.
    .. [2] Direction Choice for Accelerated Convergence in Hit-and-Run Sampling
       David E. Kaufman Robert L. Smith
       Operations Research 199846:1 , 84-95
    """
    if method == "optgp":
        sampler = GeneralizedOptGPSampler(model, processes, thinning=thinning, seed=seed)
    elif method == "achr":
        sampler = GeneralizedACHRSampler(model, thinning=thinning, seed=seed)
    else:
        raise ValueError("method must be 'optgp' or 'achr'!")

    return sampler.sample(n, fluxes = False)
