#!/bin/bash

input_file="input_CMF.txt"

# change mass in file

while read -r input_CMF; do
	echo " Processing input CMF: $input_CMF "

	sed -i '' "s/^target_CMF=[^;]*;/target_CMF=${input_CMF};/" single_model_stag_evol_CMF.m

	# change output file direcotry to follow CMF
	output_dir="CMF_sweep"

	# change these names in the file "CMF_sweep" to match the new CMF
	#sed -i '' -E "s/CMF[0-9]?_sweep/CMF_${input_CMF}_sweep/g" single_model_stag_evol_CMF.m

	# --- CMF dependent output file naming ---
	sed -i '' -E "s/(Earth_mantle_temp_CMF_)[0-9]+(\.[0-9]+)?/\1${input_CMF}/" single_model_stag_evol_CMF.m
	sed -i '' -E "s/(Earth_lith_crust_melt_CMF_)[0-9]+(\.[0-9]+)?/\1${input_CMF}/" single_model_stag_evol_CMF.m
	sed -i '' -E "s/(Earth_melt_prod_kmGyr_CMF_)[0-9]+(\.[0-9]+)?/\1${input_CMF}/" single_model_stag_evol_CMF.m
	sed -i '' -E "s/(Earth_heat_fluxWm2_CMF_)[0-9]+(\.[0-9]+)?/\1${input_CMF}/" single_model_stag_evol_CMF.m
	sed -i '' -E "s/(Earth_heat_flowTW_CMF_)[0-9]+(\.[0-9]+)?/\1${input_CMF}/" single_model_stag_evol_CMF.m

	echo ' Running model '

	arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch "single_model_stag_evol_CMF"

done < "$input_file"
