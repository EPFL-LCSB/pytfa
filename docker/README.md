# pyTFA Docker

This Docker offers a suitable environment to run pyTFA.

## Requirements

Make sure [docker](https://www.docker.com/) is installed.
You might want to install a commercial solver. Gurobi and CPLEX are supported. See [solver/instructions.txt](https://github.com/EPFL-LCSB/pytfa/blob/master/docker/solvers/instructions.txt) for more informations on how to install them.

## Running the Docker

First, build the container with `build.bat` or `. build`.
Then start the container with `run.bat` or `. run`.
```bash
. build
. run
```

You can run the examples in /pytfa/tutorials:
```bash
cd /pytfa/tutorials
python relaxation_example.py
```

You can also run them inside IPython to experiment and play with the objects:

```bash
ipython
run glycolysis_example.py
tmodel.print_info()
```
