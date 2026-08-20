#!/bin/bash

input_file_I="input_masses.txt"
input_file_II="input_radioactive_concentrations.txt"

while read -r input_mass; do
	#echo " Processing input mass: $input_mass "
	
	# --- check input mass for radius relation --
	if awk "BEGIN {exit !($input_mass > 124)}"; then
		echo " Mass > 124 M_Earth. Applying different mass-radius relation " 
		sed -i '' -E "s|^([[:space:]]*Rp=6378100\.0\*mass\^)[0-9.]+(;.*)|\10.01\2|" single_model_stag_evol.m
	else
		echo " Using standard mass scaling for radius "
		sed -i '' -E "s|^([[:space:]]*Rp=6378100\.0\*mass\^)[0-9.]+(;.*)|\10.27\2|" single_model_stag_evol.m
	fi 

	# --- mass dependent output file naming ---
	sed -i '' "s/^mass=[^;]*;/mass=${input_mass};/" single_model_stag_evol.m
	#sed -i '' -E "s/(Earth_mantle_temp_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
	#sed -i '' -E "s/(Earth_lith_crust_melt_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
	#sed -i '' -E "s/(Earth_melt_prod_kmGyr_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
	#sed -i '' -E "s/(Earth_heat_fluxWm2_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
	#sed -i '' -E "s/(Earth_heat_flowTW_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m

	while read -r input_radioactive; do
		echo " Processing input mass $input_mass at radioactive concentration: $input_radioactive"

		# --- radioactive concentration dependent output file naming ---
		sed -i '' -E "s/^U238=[0-9.]+;$/U238=${input_radioactive};/" single_model_stag_evol.m
		sed -i '' -E "s/^U235=[0-9.]+;$/U235=${input_radioactive};/" single_model_stag_evol.m
		sed -i '' -E "s/^Th=[0-9.]+;$/Th=${input_radioactive};/" single_model_stag_evol.m
		sed -i '' -E "s/^K=[0-9.]+;$/K=${input_radioactive};/" single_model_stag_evol.m
		# change output file names to include both mass and enhancement level 
		# following pattern: Earth_mantle_temp_1.00m_1.00_enhanced.dat
		sed -i '' -E "s/(Earth_mantle_temp_)[0-9.]+m(_[0-9.]+_enhanced)?\.dat/\1${input_mass}m_${input_radioactive}_enhanced.dat/" single_model_stag_evol.m
		sed -i '' -E "s/(Earth_lith_thickness_)[0-9.]+m(_[0-9.]+_enhanced)?\.dat/\1${input_mass}m_${input_radioactive}_enhanced.dat/" single_model_stag_evol.m
		sed -i '' -E "s/(Earth_lith_crust_melt_)[0-9.]+m(_[0-9.]+_enhanced)?\.dat/\1${input_mass}m_${input_radioactive}_enhanced.dat/" single_model_stag_evol.m
		sed -i '' -E "s/(Earth_melt_prod_kmGyr_)[0-9.]+m(_[0-9.]+_enhanced)?\.dat/\1${input_mass}m_${input_radioactive}_enhanced.dat/" single_model_stag_evol.m
		sed -i '' -E "s/(Earth_heat_fluxWm2_)[0-9.]+m(_[0-9.]+_enhanced)?\.dat/\1${input_mass}m_${input_radioactive}_enhanced.dat/" single_model_stag_evol.m
		sed -i '' -E "s/(Earth_heat_flowTW_)[0-9.]+m(_[0-9.]+_enhanced)?\.dat/\1${input_mass}m_${input_radioactive}_enhanced.dat/" single_model_stag_evol.m
		sed -i '' -E "s/(Earth_heat_production_)[0-9.]+m(_[0-9.]+_enhanced)?\.dat/\1${input_mass}m_${input_radioactive}_enhanced.dat/" single_model_stag_evol.m

		echo ' Running single_model_stag_evol.m '
	
		arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch "single_model_stag_evol"

	done < "$input_file_II"

done < "$input_file_I"




# reset radioactive concentrations to 1.00 m_E
echo ''
echo ' Resetting code to Earth-like parameters. '
input_radioactive=1.00
	sed -i '' -E "s/^U238=[0-9.]+;$/U238=${input_radioactive};/" single_model_stag_evol.m
	sed -i '' -E "s/^U235=[0-9.]+;$/U235=${input_radioactive};/" single_model_stag_evol.m
	sed -i '' -E "s/^Th=[0-9.]+;$/Th=${input_radioactive};/" single_model_stag_evol.m
	sed -i '' -E "s/^K=[0-9.]+;$/K=${input_radioactive};/" single_model_stag_evol.m

input_mass=1.00
sed -i '' "s/^mass=[^;]*;/mass=${input_mass};/" single_model_stag_evol.m

sed -i '' -E "s/(Earth_mantle_temp_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
sed -i '' -E "s/(Earth_lith_crust_melt_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
sed -i '' -E "s/(Earth_melt_prod_kmGyr_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
sed -i '' -E "s/(Earth_heat_fluxWm2_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
sed -i '' -E "s/(Earth_heat_flowTW_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
sed -i '' -E "s/(Earth_heat_production_)[0-9]+(\.[0-9]+)?/\1${input_radioactive}/" single_model_stag_evol.m
echo ''
echo ' End of Line '
