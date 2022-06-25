#!/bin/bash

input_file=$1

make

while read line;
do

	echo $line
	./main.x $line 

done < <(grep -v '#' $input_file)

rm fit_square.txt
mv bound*.txt DATA/BOUND/
