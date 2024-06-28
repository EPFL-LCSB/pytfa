#!/bin/sh
#

echo "Installing Gurobi"

if [ -f /solvers/gurobi9.0.1_linux64.tar.gz ]; then \
    cp /solvers/gurobi9.0.1_linux64.tar.gz /opt &&  \
    cd /opt                     &&  \
    tar xvfz gurobi9.0.1_linux64.tar.gz     &&  \
    cd /opt/gurobi901/linux64           &&  \
    python3 setup.py install            &&  \
    rm /opt/gurobi9.0.1_linux64.tar.gz; fi

# if [ -f /solvers/gurobi8.1.1_linux64.tar.gz ]; then \
#     cp /solvers/gurobi8.1.1_linux64.tar.gz /opt &&  \
#     cd /opt                     &&  \
#     tar xvfz gurobi8.1.1_linux64.tar.gz     &&  \
#     cd /opt/gurobi811/linux64           &&  \
#     python3 setup.py install            &&  \
#     rm /opt/gurobi8.1.1_linux64.tar.gz; fi

#if [ -f /solvers/gurobi8.0.0_linux64.tar.gz ]; then \
#	cp /solvers/gurobi8.0.0_linux64.tar.gz /opt	&&	\
#	cd /opt						&&	\
#	tar xvfz gurobi8.0.0_linux64.tar.gz		&&	\
#	cd /opt/gurobi800/linux64 			&& 	\
#	python3 setup.py install			&&	\
#	rm /opt/gurobi8.0.0_linux64.tar.gz; fi

#if [ -f /solvers/gurobi7.5.2_linux64.tar.gz ]; then \
#	cp /solvers/gurobi7.5.2_linux64.tar.gz /opt	&&	\
#	cd /opt						&&	\
#	tar xvfz gurobi7.5.2_linux64.tar.gz		&&	\
#	cd /opt/gurobi752/linux64 			&& 	\
#	python3 setup.py install			&&	\
#	rm /opt/gurobi7.5.2_linux64.tar.gz; fi
	

#if [ -f /solvers/gurobi7.0.2_linux64.tar.gz ]; then \
#	cp /solvers/gurobi7.0.2_linux64.tar.gz /opt	&&	\
#	cd /opt						&&	\
#	tar xvfz gurobi7.0.2_linux64.tar.gz		&&	\
#	cd /opt/gurobi702/linux64 			&& 	\
#	python3 setup.py install			&&	\
#	rm /opt/gurobi7.0.2_linux64.tar.gz; fi
#	
