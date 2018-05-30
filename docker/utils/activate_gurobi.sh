#!/bin/sh
# Choose the GHOME variable depending on your gurobi version

export GHOME_702="/opt/gurobi702/linux64"
export GHOME_752="/opt/gurobi752/linux64"
export GHOME_800="/opt/gurobi800/linux64"
export GHOME=$GHOME_800

if [ -d $GHOME ]; then \
	#$GHOME/bin/grbgetkey $LICENSE_KEY	&& \

	# Export path
	echo "export GUROBI_HOME=${GHOME}" >> ~/.bashrc			&&	\
	echo 'export PATH="${PATH}:${GUROBI_HOME}/bin"' >> ~/.bashrc			&&	\
	echo 'export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"' >> ~/.bashrc; fi
