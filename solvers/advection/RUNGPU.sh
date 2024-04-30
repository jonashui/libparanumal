#!/bin/sh
# embedded options to bsub - start with #BSUB
# -- our name ---
#BSUB -J SWEtest
# -- choose queue --
#BSUB -q gpuv100
# -- Notify me by email when execution begins --
#BSUB -B
# -- Notify me by email when execution ends   --
#BSUB -N
# -- email address -- 
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
##BSUB -u your_email_address
# -- Output File --
#BSUB -o Output_%J.txt
# -- Error File --
#BSUB -e Error_%J.txt

#BSUB -n 4
#BSUB -R "rusage[mem=15GB]"
#BSUB -R "span[ptile=1]"
#BSUB -W 1:00
#BSUB -gpu "num=1"

# Load modules
module purge
module load mpi/3.1.3-gcc-8.2.0
module load gcc/12.2.0-binutils-2.39
module load openblas
module load cuda/10.0

# Needed environment variables
export OCCA_DIR=~/libparanumal/occa
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

# Make occa for CPU - Only needed if already compiled for GPU or not compiled yet
cd ../../occa
make clean
make -j
cd ../solvers/advection/

# Build project
#make clean
make -j `nproc` 
#./advectionMain setups/setupTri2D.rc
mpirun ./advectionMain setups/setupTri2D.rc