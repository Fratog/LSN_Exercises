#!/bin/bash

input_file=$1

make

#https://stackoverflow.com/questions/43741612/mpirun-breaks-out-the-while-loop-in-bash
#https://it.wikipedia.org/wiki//dev/null
while read line
	do

		echo $line
		mpiexec -np 8 main.out $line > /dev/null < /dev/null

	done < <(grep -v '#' $input_file)



