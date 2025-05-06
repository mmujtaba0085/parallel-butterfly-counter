# Performance Analysis Guide

Note: use ssh -X mpi@main (to login, instead of su - mpi), this is to also get GUI analysis

### TAU SETUP

```bash
sudo apt update
sudo apt install build-essential gfortran g++ make cmake wget git libpapi-dev libotf-dev libelf-dev


cd ~
git clone https://github.com/UO-OACISS/tau2.git
cd tau2
./configure -mpi -openmp

```

### Build TAU

```bash
make -j$(nproc)
```

At this stage, we'll find tau_cxx.sh and paraprof inside tau2/x86_64/bin/

### Add/Export ENVs

this is for one time setup (permanent), otherwises you'll have to run these commands everytime, thus we'll make it easier by exporting and adding to bash profiles

```bash
export TAU_MAKEFILE=$HOME/tau2/x86_64/lib/Makefile.tau-mpi-openmp
export PATH=$HOME/tau2/x86_64/bin:$PATH

echo 'export TAU_MAKEFILE=$HOME/tau2/x86_64/lib/Makefile.tau-mpi-openmp' >> ~/.bashrc
echo 'export TAU_MAKEFILE=$HOME/tau2/x86_64/lib/Makefile.tau-mpi-openmp' >> ~/.bash_profile
echo 'export TAU_MAKEFILE=$HOME/tau2/x86_64/lib/Makefile.tau-mpi-openmp' >> ~/.profile
echo 'export PATH=$HOME/tau2/x86_64/bin:$PATH' >> ~/.bashrc
echo 'export PATH=$HOME/tau2/x86_64/bin:$PATH' >> ~/.bash_profile
echo 'export PATH=$HOME/tau2/x86_64/bin:$PATH' >> ~/.profile

source ~/.bashrc
source ~/.bash_profile
source ~/.profile
```

### Compile and Run Command (From Root)

```bash
tau_cxx.sh -openmp -O2 parallel-butterfly-counter/Parallel/MPI+OpenMP.cpp -o parallel-butterfly-counter/hybrid_butterfly -lmetis //compile

mpiexec -n 4 -hosts master,slave ./parallel-butterfly-counter/hybrid_butterfly parallel-butterfly-counter/datasets/bipartite_graph_100k.txt

```

after running the program, profiles will be generated in the root folder (ls to check)
just run the following command to view these

```bash
pprof
paraprof
```

### pprof -> text/tabular_analysis

### paraprof -> gui based

---

## Better Performance Analysis with Automation

### analyze.sh

```bash
#!/bin/bash

# Define directories and filenames
SOURCE_FILE="parallel-butterfly-counter/Parallel/MPI+OpenMP.cpp"
OUTPUT_BINARY="parallel-butterfly-counter/hybrid_butterfly"
DATA_FILE="datasets/bipartite_graph_100k.txt"

# compile the code with TAU
echo "Compiling code with TAU..."
tau_cxx.sh -openmp -O2 $SOURCE_FILE -o $OUTPUT_BINARY -lmetis

# run the program using MPI
echo "Running the program using MPI..."
mpiexec -n 4 -hosts master,slave ./$OUTPUT_BINARY $DATA_FILE

# performance analysis with pprof
echo "Showing performance analysis..."
pprof

# performance analysis with paraprof
echo "Showing performance analysis..."
paraprof

# Uncomment the line below if you want to filter for MPI overhead
# pprof -s | grep MPI

echo "Performance analysis complete."

```
