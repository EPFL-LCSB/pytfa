# pyTFA Docker

This Docker offers a suitable environment to run Thermodynamic constraints for flux-based reactions, in Python.

## Requirements

This module requires [EPFL-LCSB/pytfa](https://github.com/EPFL-LCSB/pytfa) to work properly. You can activate it by running (inside the pytfa folder repository)
```bash
git submodule init
git submodule update
```
You might also want to install a dedicated solver, GLPK or CPLEX for example.

## Running the Docker

Make sure [docker](https://www.docker.com/) is installed, then build the container first
with `build.bat` or `build.sh`.

Then start the container with `run.bat` or `run.sh`.

