#!/bin/bash

phase=$1

cp "input."$phase".eq" "input.in"
./MD.x

mv *.dat "DATA/"$phase"_equilibration/"
cp *.out "DATA/"$phase"_equilibration/"


cp "input."$phase "input.in"
./MD.x

mv *.dat "DATA/"$phase"_averages/"
cp *.out "DATA/"$phase"_averages/"

rm *.out
