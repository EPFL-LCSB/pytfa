#!/bin/sh
# Thanks to cdiener: https://hub.docker.com/r/cdiener/cobra-docker/~/dockerfile/
# For the solution of simply getting the bins and python hooks

#

echo "Installing and Moving CPLEX files"

# Default Py3.5 install
#if [ -d /solvers/ibm ]; then cd /solvers/ibm/ILOG/CPLEX_Studio1271/cplex/python/3.5/x86-64_linux/ && \#
#	python3 setup.py install && \
#	cp /solvers/ibm/ILOG/CPLEX_Studio1271/cplex/bin/x86-64_linux/cplex /usr/bin/; fi

# Default Py3.6 install
if [ -d /solvers/ibm ]; then cd /solvers/ibm/ILOG/CPLEX_Studio128/cplex/python/3.6/x86-64_linux/ && \
	python3 setup.py install && \
	cp /solvers/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/cplex /usr/bin/; fi

# Tweaked Py2.6/2.7 install on conda Venv

# if [ -d /solvers/ibm ]; then cd /solvers/ibm/cplex/python/2.6/x86-64_linux/ && \
	# python setup.py install && \
	# cp /solvers/ibm/cplex/bin/x86-64_linux/cplex /usr/bin/; fi
