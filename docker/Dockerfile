FROM python:3.6
# Warning: cplex 12.7.1 is only compatible with python 3.5, while
# Gurobi 7.5 is only compatbile with with python 3.6. If you wan to
# install both, I recommend downgrading the gurobi version.
# Another way (not recommended) is tweaking the setup.py python version
# requirements of either one of them, but this might alter stability
# of the solvers

# Install missing deps
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2-dev     \
        libxslt1-dev    \
        less            \
    && rm -rf /var/lib/apt/lists/*

ENV USER pytfa
ENV HOME /home/$USER

RUN useradd -ms "/bin/bash" "$USER"
USER $USER
WORKDIR $HOME

USER root

# Copy python package requirements
COPY requirements.txt .

# Install python packages
RUN pip install -r requirements.txt

# Take care of the solvers
COPY ./solvers /solvers
COPY ./utils /utils

RUN chmod u+x /utils/*.sh

# Install CPLEX
RUN /utils/install_cplex.sh
# Install gurobi
COPY ./utils/gurobi.lic* ./
RUN /utils/install_gurobi.sh

# Remove installers
RUN rm -rf /solvers

# Add extra src if needed
RUN mkdir /src
COPY src/ /src/

# Make the /src/pytfa folder that will link the sources
RUN mkdir /src/pytfa

COPY .bashrc $HOME
RUN chown "$USER" "$HOME/.bashrc"

#Finalizing installation

RUN chmod +x /utils/activate_gurobi.sh 

USER $USER
RUN mkdir ./work

# Activation of necessary licenses
RUN /utils/activate_gurobi.sh

# Load your package in development mode on startup
ENTRYPOINT ["/bin/bash", "-c", "pip install --user -e /src/pytfa && $0 $*"]
CMD /bin/bash
