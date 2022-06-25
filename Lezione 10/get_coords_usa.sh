#!/bin/bash

coords="American_capitals.dat" #American_capitals.dat
output="usa_only_coords.dat"

#create a suitable file for capitals' coordinates
awk '{ print $(NF-1), $NF}' $coords > $output 
tail -n +2  $output > "appo.tmp" && mv "appo.tmp" $output
