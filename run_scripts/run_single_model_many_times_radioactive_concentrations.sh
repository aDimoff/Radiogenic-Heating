#!/bin/bash

input_file="input_radioactive_concentrations.txt"
	
while read -r input_radioactive; do
	echo " Processing input radioactive concentration: $input_radioactive "
	# --- radioactive concentration dependent output file naming ---
	#sed -i '' -E "s/^U238=[0-9.]+;$/U238=${input_radioactive};/" single_model_stag_evol.m
	#sed -i '' -E "s/^U235=[0-9.]+;$/U235=${input_radioactive};/" single_model_stag_evol.m
	#sed -i '' -E "s/^Th=[0-9.]+;$/Th=${input_radioactive};/" single_model_stag_evol.m
	sed -i '' -E "s/^K=[0-9.]+;$/K=${input_radioactive};/" single_model_stag_evol.m
	sed -i '' -E "s/(Earth_mantle_temp_K_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
	sed -i '' -E "s/(Earth_lith_thickness_K_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
	sed -i '' -E "s/(Earth_lith_crust_melt_K_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
	sed -i '' -E "s/(Earth_melt_prod_kmGyr_K_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
	sed -i '' -E "s/(Earth_heat_fluxWm2_K_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
	sed -i '' -E "s/(Earth_heat_flowTW_K_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
	echo ' Running single_model_stag_evol.m '

	#arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch "single_model_stag_evol"
done < "$input_file"

# reset radioactive concentrations to 1.000 m_E
echo ''
echo ' Resetting code to Earth-like parameters. '
input_radioactive=1.000
	sed -i '' -E "s/^U238=[0-9.]+;$/U238=${input_radioactive};/" single_model_stag_evol.m
	sed -i '' -E "s/^U235=[0-9.]+;$/U235=${input_radioactive};/" single_model_stag_evol.m
	sed -i '' -E "s/^Th=[0-9.]+;$/Th=${input_radioactive};/" single_model_stag_evol.m
	sed -i '' -E "s/^K=[0-9.]+;$/K=${input_radioactive};/" single_model_stag_evol.m

echo ''
echo ' End of Line '
