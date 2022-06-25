#!/bin/bash

#this script must ensure that all the equilbration part is done, and there are the ouputfiles necessary to execute 

#get input values
Type_of_ev=$1  #MC for monte carlo , MD for dynamic evolution
Nstep=$2 # will be the number of istanenteuneus values to use for equilibratoin, energy autocorellation data...
restart=0 #1 if true 0 if false


make

phases=("solid" "liquid" "gas") 
for phase in ${phases[@]}; 
do

	if [ "$Type_of_ev" == "MC" ]
	then
		bool=1 #flag to pass to the main
		path=DATA/$phase
	else
		bool=0
		path=DATA/Molecular_Dyn
	fi
echo $bool

	#Nstep is actually the number of blocks nblk and the 1 after is nstep. This way we can get istantaneous values.
	cat "input."$phase | awk -v tofev=$bool -v Nstep=$Nstep -v restart=$restart 'FNR > 1 {print(tofev, restart, $3, $4, $5, $6, $7, Nstep, 1, $10)}' > "input.in"
	./Molecular.x #the executable takes an input which is the number of moves to equilibrate the system immediately after initalization. Here we study the equilibratoin process, so for now it is 0.
	
	rm gr_final_*
	rm *.out
	mv output_*.dat DATA/Exercise7.2/Equilibration/

# you can use restart bool to decide if move files in Equilivration path or in Exercise path, do it, and make exercie7.2 a variabile	
	#mv output* DATA/Exercise7.3/
	#mv gr_final* DATA/Exercise7.3/
	#cp config* DATA/Exercise7.3/
	#change names of equilbirum files

#move the various aoutputs in the correct way
done 
 
#########file paths




