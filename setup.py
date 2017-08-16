""" Thermodynamic constraints for Flux-Based Analysis of reactions

.. moduleauthor:: pyTFA team


"""

from setuptools import setup

setup(name='pytfa',
      version='1.1.0',
      author='pyTFA team',
      install_requires=['cobra>0.6,<=0.8.1','bokeh>=0.12.1','optlang'],
      py_modules=['pytfa'],
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='pyTFA adds thermodynamic constraints for Flux-Based Analysis of reactions'
     )
