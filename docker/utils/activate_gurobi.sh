#!/bin/sh
# Script to activate gurobi - remember to add your license key!\\
# Rename to activate_gurobi.sh
# Gurobi does not allow to install a single-user license on root
# Choose the GHOME variable depending on your gurobi version

export LICENSE_KEY=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
export GHOME_702="/opt/gurobi702/linux64"
export GHOME_750="/opt/gurobi750/linux64"
export GHOME=$GHOME_702

if [ -d $GHOME ]; then \
	$GHOME/bin/grbgetkey $LICENSE_KEY	&& \

	# Export path
	echo "export GUROBI_HOME=${GHOME}" >> ~/.bashrc			&&	\
	echo 'export PATH="${PATH}:${GUROBI_HOME}/bin"' >> ~/.bashrc			&&	\
	echo 'export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"' >> ~/.bashrc; fi
