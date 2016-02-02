#!/bin/bash
module purge
module load null
export PATH=/short/v29/mma659/bin/:/apps/gcc/4.6.4/bin/:$HOME/bin/:$PATH
echo PATH=$PATH
export LD_LIBRARY_PATH=/short/v29/mma659/lib/:/apps/fftw3/3.3.3/lib/:/apps/gcc/4.6.4/lib64/:$HOME/lib/:/usr/lib/:/usr/lib64/:$LD_LIBRARY_PATH
echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export MACHINE=raijin_cluster
echo MACHINE=$MACHINE
