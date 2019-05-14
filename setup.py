""" Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team


"""

from setuptools import setup, find_packages
# import os
# from pip.req import parse_requirements
# from pip.download import PipSession

# __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
#
# def read_requirements():
#     '''parses requirements from requirements.txt'''
#     reqs_path = os.path.join(__location__, 'requirements.txt')
#     install_reqs = parse_requirements(reqs_path, session=PipSession())
#     reqs = [str(ir.req) for ir in install_reqs]
#     return reqs


version_tag = '0.9.0-b1'

setup(name='pytfa',
      version=version_tag,
      author='pyTFA team',
      author_email='softwares.lcsb@epfl.ch',
      url='https://github.com/EPFL-LCSB/pytfa/',
      download_url='https://github.com/EPFL-LCSB/pytfa/archive/'+version_tag+'.tar.gz',
      install_requires=['cobra>0.11',
                        'bokeh>=0.12.1',
                        'networkx',
                        'optlang',
                        'pytest',
                        'scipy',
                        'tqdm'],
      packages = find_packages(),
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='pyTFA adds Thermodynamics-based Flux Analysis',
      keywords=['pytfa','tfa','thermodynamics','flux analysis'],

      license='Apache 2.0',

      # See https://PyPI.python.org/PyPI?%3Aaction=list_classifiers
      classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 4 - Beta',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Environment :: Console',

            # Pick your license as you wish (should match "license" above)
            'License :: OSI Approved :: Apache Software License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
      ],
     )
