Solver setup
============

This document is written assuming a Docker container installation. 
However, you can easily adapt the content to other types of Linux-based installations.

GPLK
-----

GLPK should be directly available from the requirements.

CPLEX
-----

You will need to first install CPLEX on a Linux machine.

Place in `etfl/docker/solvers/` the folder `/ibm` that is installed by CPLEX (usually in `/opt`).
You actually only need the following substructure (makes the container lighter):

.. code-block:: text

    .
    └───ibm
        └───ILOG
            └───CPLEX_StudioXXXX
                └───cplex
                    ├───bin
                    ├───include
                    ├───lib
                    └───python
				
Gurobi
------

Place in `etfl/docker/solvers/` the tarball you downloaded from the website, and modify accordingly the files:

.. code-block:: text

	../utils/install_gurobi.sh 
	../utils/activate_gurobi.sh

Make sure you change the paths and filenames to reflect the actual version of Gurobi you are running.

Gurobi needs a floating license for Docker instances, (see http://www.gurobi.com/documentation/7.5/quickstart_windows/setting_up_and_using_a_flo.html#subsection:tokenserver)
Once your system administrator set it up, you will need to add your gurobi license server to ../utils/gurobi.lic.template, and rename it to gurobi.lic
