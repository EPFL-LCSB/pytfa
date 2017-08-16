# pyTFA Docker

This Docker offers a suitable environment to run Thermodynamic constraints for flux-based reactions, in Python.

## Requirements

Make sure [docker](https://www.docker.com/) is installed.
You might want to install a dedicated solver, GLPK or CPLEX for example.

## Running the Docker

First, build the container with `build.bat` or `build.sh`.
Then start the container with `run.bat` or `run.sh`.

You can run the examples in /pytfa/tutorials:
```bash
cd /pytfa/tutorials
python glycolysis_example.py
```

You can also run them inside IPython to experiment:

```bash
ipython
run glycolysis_example.py
tmodel.print_info()
```