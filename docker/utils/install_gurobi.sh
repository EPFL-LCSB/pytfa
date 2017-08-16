#!/bin/sh
#

echo "Installing Gurobi"

# if [ -f /solvers/gurobi7.5.0_linux64.tar.gz ]; then \
# 	cp /solvers/gurobi7.5.0_linux64.tar.gz /opt	&&	\
# 	cd /opt						&&	\
# 	tar xvfz gurobi7.5.0_linux64.tar.gz		&&	\
# 	cd /opt/gurobi750/linux64 			&& 	\
# 	python3 setup.py install			&&	\
# 	rm gurobi7.0.2_linux64.tar.gz; fi
	

if [ -f /solvers/gurobi7.0.2_linux64.tar.gz ]; then \
	cp /solvers/gurobi7.0.2_linux64.tar.gz /opt	&&	\
	cd /opt						&&	\
	tar xvfz gurobi7.0.2_linux64.tar.gz		&&	\
	cd /opt/gurobi702/linux64 			&& 	\
	python3 setup.py install			&&	\
	rm /opt/gurobi7.0.2_linux64.tar.gz; fi
	
