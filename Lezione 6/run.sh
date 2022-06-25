#!/bin/bash

input=$1 #input file
samp=$2 #M for metropolis, G for gibbs
h=$3 #magnetic field


make
#########file paths
if [ "$samp" == "M" ]
then
	bool=1
	path=DATA/Metropolis
else
	bool=0
	path=DATA/Gibbs
fi

##########create input file
cat $input | awk -v s=$bool -v f=$h 'FNR > 1 {print($1, $2, $3, f, s, $6, $7, $8)}' > "my_"$input


##########run simulations
input="my_"$input
while read line;
do

	echo $line
	./Monte_Carlo_ISING_1D.x $line

done < <(grep -v '#' $input)

########movr files and get final results
	mv $input $path
	#mv "equil"*$samp".txt" $eq_path
	mv *$samp".txt" $path
	rm config* ##just for now

#https://stackoverflow.com/questions/8654051/how-can-i-compare-two-floating-point-numbers-in-bash
if [ $h == 0 ] #using == for numbers is technincally wrong, but bash doesn't work floating points. with this implementation we avoid problems. BUT if you want to work at zero external field we MUST pass 0 to h and NOT 0.0 otherwise string comparison fails. MODIFY IT
then
	echo "here"	
	obs=("Energy" "Heat" "Chi")
else
	obs=("Mag")
fi

for str in ${obs[@]}; 
do

	for i in $(ls $path/$str*.txt)
	do
		cat $i | tail -n 1 >> $str"_results_"$samp".txt"
	done
done

mv *$samp".txt" $path

