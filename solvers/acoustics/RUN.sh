#!/bin/sh
# embedded options to bsub - start with #BSUB
# -- our name ---
#BSUB -J SWEtest
# -- choose queue --
#BSUB -q hpc
# -- Notify me by email when execution begins --
#BSUB -B
# -- Notify me by email when execution ends   --
#BSUB -N
# -- email address -- 
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
##BSUB -u your_email_address
# -- Output File --
#BSUB -o aOutput_%J.txt
# -- Error File --
#BSUB -e aError_%J.txt

#BSUB -n 4
#BSUB -R "rusage[mem=10GB]"
#BSUB -R "span[block=1]"
#BSUB -W 1:00
#BSUB -R "select[model == XeonGold6142]"

# Load modules
module purge
module load mpi/3.1.3-gcc-8.2.0
module load gcc/12.2.0-binutils-2.39
module load openblas

# Needed environment variables
#export OCCA_DIR=~/libparanumal/occa
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

# Make occa for CPU - Only needed if already compiled for GPU or not compiled yet
#cd ../../occa
#make clean
#make -j
#cd ../solvers/advection/

# Build project
#make realclean
make clean
make -j
#./advectionMain setups/setupTri2D.rc
mpiexec -np 1 --map-by slot:PE=4 ./acousticsMain setups/setupTri2D.rc