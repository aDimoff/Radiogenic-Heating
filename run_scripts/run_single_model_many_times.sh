#!/bin/bash

input_file="input_masses.txt"

while read -r input_mass; do
	echo " Processing input mass: $input_mass "
	
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
	sed -i '' -E "s/(Earth_mantle_temp_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
	sed -i '' -E "s/(Earth_lith_thickness_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
	sed -i '' -E "s/(Earth_lith_crust_melt_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
	sed -i '' -E "s/(Earth_melt_prod_kmGyr_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
	sed -i '' -E "s/(Earth_heat_fluxWm2_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
	sed -i '' -E "s/(Earth_heat_flowTW_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m

	echo ' Running single_model_stag_evol.m '
	
	arch -x86_64 /Applications/MATLAB_R2020a.app/bin/matlab -batch "single_model_stag_evol"

done < "$input_file"


# reset mass to 1.00 m_E
echo ''
echo ' Resetting code to Earth-like parameters. '
input_mass=1.00
sed -i '' "s/^mass=[^;]*;/mass=${input_mass};/" single_model_stag_evol.m
sed -i '' -E "s/(Earth_mantle_temp_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
sed -i '' -E "s/(Earth_lith_thickness_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
sed -i '' -E "s/(Earth_lith_crust_melt_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
sed -i '' -E "s/(Earth_melt_prod_kmGyr_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
sed -i '' -E "s/(Earth_heat_fluxWm2_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m
sed -i '' -E "s/(Earth_heat_flowTW_)[0-9]+(\.[0-9]+)?m/\1${input_mass}m/" single_model_stag_evol.m

# reset mass - radius relation to original
sed -i '' -E "s|^([[:space:]]*Rp=6378100\.0\*mass\^)[0-9.]+(;.*)|\10.27\2|" single_model_stag_evol.m

echo ''
echo ' End of Line '
