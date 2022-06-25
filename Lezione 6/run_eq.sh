#!/bin/bash

nstep=$1 #it is actually nblk and nstep=1. this way we get istantenous values

./Monte_Carlo_ISING_1D.x 0.5 50 1.0 0.0 1 $nstep 1 0
mv Energy*.txt DATA/Eq_Metropolis
rm *.txt

./Monte_Carlo_ISING_1D.x 0.5 50 1.0 0.0 0 $nstep 1 0
mv Energy*.txt DATA/Eq_Gibbs
rm *.txt

#rm config*
