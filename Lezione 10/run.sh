#!/bin/bash

#compile
make

flag=$1 #must be 1
row1=$2 # here row 1 and row2 can be const, since row1 const means we work with fixed coords in square (here is ok, we are not stimating E[dmin]
row2=$3 #even if row2 is const in the main we have row2+rank, so the differnt cores will have anyway a statistically indpendent evolution
Ncities=$4 
selection_index=$5 #2 for standard 1 for SA selection operator

mpiexec -np 1 main.out 0 $row1 $row2 $Ncities $selection_index #no parallelization

for i in {2..8}
	do
		mpiexec -np $i main.out $flag $row1 $row2 $Ncities $selection_index #parallelization with varying number of cores
	done
	
if [ $5 -eq 1 ]
	then
		outfile="fit_varying_ncores_"$Ncities"_SA_nonpar.txt"
elif [ $5 -eq 2 ]
	then
		outfile="fit_varying_ncores_"$Ncities"_Standard_nonpar.txt"
else
	echo "ERROR in run.sh"
	fi

#echo $outfile

if [ $flag -eq 1 ]
	then	
	./kfact.sh $outfile $Ncities	
	fi	
