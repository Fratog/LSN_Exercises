#!/bin/bash

input=$1
Ncities=$2

output="kfact_"$Ncities".txt"

awk -v fname=$output '{if(NR==1) appo=$2;}{print($1, (appo-$2)/$2) >> fname}' $input

mv $input DATA/
mv $output DATA/
