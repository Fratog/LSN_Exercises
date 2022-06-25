#!/bin/bash

#get input values
Type_of_ev=$1  #MC for monte carlo , MD for dynamic evolution
Nblocks=$2
Nstep=$3 
restart=0 # i'll put the equilibration in the main, so restar won't be necessary, put it to 0


make 

phases=("solid" "liquid" "gas") 
for phase in ${phases[@]}; 
do

	if [ "$Type_of_ev" == "MC" ]
	then
		bool=1 #flag for Monte carlo or Verlet option to pass to the input file
		input="input."$phase 
	else
		bool=0
		input="input."$phase".verlet"
	fi

	#modify input according to values passed to the script
	cat "$input" | awk -v tofev=$bool -v Nblocks=$Nblocks -v Nstep=$Nstep -v restart=$restart 'FNR > 1 {print(tofev, restart, $3, $4, $5, $6, $7, Nblocks, Nstep, $10)}' > "input.in"
	./Molecular.x 

#move the various outputs in the coorect way	
done 

