#!/bin/sh
# embedded options to bsub - start with #BSUB
# -- our name ---
#BSUB -J SWEtest
# -- choose queue --
#BSUB -q gpua100
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

#BSUB -gpu "num=1:mode=exclusive_process"
##BSUB -gpu "num=1"

#BSUB -n 4
#BSUB -R "rusage[mem=5GB]"
#BSUB -R "span[block=1]"
#BSUB -W 15:00

# Load modules
module purge
module load mpi/4.0.5-gcc-8.4.0
module load openblas
module load cuda/12.0

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
#make clean
#make -j
#./advectionMain setups/setupTri2D.rc
mpiexec -np 1 --map-by slot:PE=4 ./SWECMain setups/setupTri2D.rc