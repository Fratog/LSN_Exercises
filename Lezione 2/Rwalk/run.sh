#!/bin/bash

make

while IFS=" " read -r Nstep Nsim Nblock dim
do

./RW.x $Nstep $Nsim $Nblock $dim

done < input

mv *.txt DATA/

