Integrating metabolomics
========================

In this short example we will go through a simple case of integration of absolute metabolite concentrations.

Let us imagine we got absolute concentration values for cytosolic ATP:

.. math::

	5 . 10^{-3} \; mol.L^{-1} \le [X] \le 3 . 10^{-2} \; mol.L^{-1}


Then you can tell the model that your (log) concentration is limited in range:

.. code:: python

	from math import log

	mymodel.log_concentration.atp_c.variable.lb = log(5e-3)
	mymodel.log_concentration.atp_c.variable.ub = log(3e-2)


This will constrain the dG according to your concentration measurements for cytosolic ATP.
As a reminder, the dG (not the dGo) takes activity (here, concentrations) into account for its calculation. You can find a more detailed explanation in those papers:

* Henry, Christopher S., Linda J. Broadbelt, and Vassily Hatzimanikatis. "Thermodynamics-based metabolic flux analysis." Biophysical journal 92.5 (2007): 1792-1805.
    
* Soh, Keng Cher, Ljubisa Miskovic, and Vassily Hatzimanikatis. "From network models to network responses: integration of thermodynamic and kinetic properties of yeast genome-scale metabolic networks." FEMS yeast research 12.2 (2012): 129-143.

