""" Thermodynamic constraints for Flux-Based Analysis of reactions

.. moduleauthor:: pyTFA team


"""

from setuptools import setup

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='pytfa',
      version='1.1.0',
      author='pyTFA team',
      install_requires=requirements,
      py_modules=['pytfa'],
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='pyTFA adds thermodynamic constraints for Flux-Based Analysis of reactions'
     )
