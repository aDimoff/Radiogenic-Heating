#!/bin/bash

input_file="input_masses.txt"

# change mass in file

while read -r input_mass; do
	echo " Processing input mass: $input_mass "

	sed -i '' "s/^mass=[^;]*;/mass=${input_mass};/" single_model_stag_evol.m

	# change output file direcotry to follow mass
	output_dir="mass_${input_mass}_enhanced_sweep"

	# change these names in the file "mass4_enhanced_sweep" to match the new mass
	# making sure the input mass is an integer for the file name
	sed -i '' -E "s/mass[0-9]?_enhanced_sweep/mass_${input_mass}_enhanced_sweep/g" single_model_stag_evol.m

	echo ' Running sweep '

	./run_single_model_many_times_radioactive_concentrations.sh

done < "$input_file"
