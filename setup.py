""" Thermodynamic constraints for Flux-Based Analysis of reactions

.. moduleauthor:: pyTFA team


"""

from setuptools import setup
from pip.req import parse_requirements
from pip.download import PipSession

def read_requirements():
    '''parses requirements from requirements.txt'''
    reqs_path = os.path.join(__location__, 'requirements.txt')
    install_reqs = parse_requirements(reqs_path, session=PipSession())
    reqs = [str(ir.req) for ir in install_reqs]
    return reqs

setup(name='pytfa',
      version='0.6.0',
      author='pyTFA team',
      author_email='softwares.lcsb@epfl.ch',
      url='https://github.com/EPFL-LCSB/pytfa/',
      download_url='https://github.com/EPFL-LCSB/pytfa/archive/0.6.0-a1.tar.gz',
      install_requires=read_requirements(),
      py_modules=['pytfa'],
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='pyTFA adds thermodynamic constraints for Flux-Based Analysis of reactions',
      keywords=['pytfa','tfa','thermodynamics','flux analysis'],
     )
