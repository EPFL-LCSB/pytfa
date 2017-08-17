""" Thermodynamic constraints for Flux-Based Analysis of reactions

.. moduleauthor:: pyTFA team


"""

from setuptools import setup
import os
from pip.req import parse_requirements
from pip.download import PipSession

# __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
#
# def read_requirements():
#     '''parses requirements from requirements.txt'''
#     reqs_path = os.path.join(__location__, 'requirements.txt')
#     install_reqs = parse_requirements(reqs_path, session=PipSession())
#     reqs = [str(ir.req) for ir in install_reqs]
#     return reqs

version_tag = '0.6.1-a3'

setup(name='pytfa',
      version=version_tag,
      author='pyTFA team',
      author_email='softwares.lcsb@epfl.ch',
      url='https://github.com/EPFL-LCSB/pytfa/',
      download_url='https://github.com/EPFL-LCSB/pytfa/archive/'+version_tag+'.tar.gz',
      install_requires=['cobra>0.6',
                        'bokeh>=0.12.1',
                        'optlang',
                        'pytest',
                        'scipy'],
      py_modules=['pytfa'],
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='pyTFA adds thermodynamic constraints for Flux-Based Analysis of reactions',
      keywords=['pytfa','tfa','thermodynamics','flux analysis'],
     )
