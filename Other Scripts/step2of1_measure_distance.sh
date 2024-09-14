#!/bin/bash

# Define the pairs of residues and atoms
declare -A pairs
# Example 
# pairs["D69_OD1_N300_ND2"]="resid 69 and name OD1 plus resid 300 and name ND2"
# pairs [""]=

# Layer1
pairs["F260_CA_I111_CA"]="resid 260 and name CA plus resid 111 and name CA"
pairs["W264_CA_I111_CA"]="resid 264 and name CA plus resid 111 and name CA"
pairs["G203_CA_F260_CA"]="resid 203 and name CA plus resid 260 and name CA"

# Layer2
pairs["V256_CA_L114_CA"]="resid 256 and name CA plus resid 114 and name CA"
pairs["M257_CA_L114_CA"]="resid 257 and name CA plus resid 114 and name CA"
pairs["V256_CA_Y210_CA"]="resid 256 and name CA plus resid 210 and name CA"
pairs["M257_CA_L207_CA"]="resid 257 and name CA plus resid 207 and name CA"
pairs["N300_ND2_L114_CA"]="resid 300 and name ND2 plus resid 114 and name CA"

# Layer3
pairs["Y304_CA_V44_CA"]="resid 304 and name CA plus resid 44 and name CA"
pairs["Y304_CA_F311_CA"]="resid 304 and name CA plus resid 311 and name CA"
pairs["Y304_CA_L114_CA"]="resid 304 and name CA plus resid 114 and name CA"
pairs["Y304_CA_I117_CA"]="resid 304 and name CA plus resid 117 and name CA"
pairs["Y304_CA_R121_NH1"]="resid 304 and name CA plus resid 121 and name NH1"
pairs["Y304_CA_R121_NH2"]="resid 304 and name CA plus resid 121 and name NH2"

# Layer4
pairs["R121_NH1_D120_OD1"]="resid 121 and name NH1 plus resid 120 and name OD1"
pairs["R121_NH1_D120_OD2"]="resid 121 and name NH1 plus resid 120 and name OD2"
pairs["R121_NH2_D120_OD1"]="resid 121 and name NH2 plus resid 120 and name OD1"
pairs["R121_NH2_D120_OD2"]="resid 121 and name NH2 plus resid 120 and name OD2"
pairs["R121_NH1_E246_OE1"]="resid 121 and name NH1 plus resid 246 and name OE1"
pairs["R121_NH1_E246_OE2"]="resid 121 and name NH1 plus resid 246 and name OE2"
pairs["R121_NH2_E246_OE1"]="resid 121 and name NH2 plus resid 246 and name OE1"
pairs["R121_NH2_E246_OE2"]="resid 121 and name NH2 plus resid 246 and name OE2"
pairs["R121_NH1_L253_CA"]="resid 121 and name NH1 plus resid 253 and name CA"
pairs["R121_NH2_L253_CA"]="resid 121 and name NH2 plus resid 253 and name CA"
pairs["R121_NH1_A124_CA"]="resid 121 and name NH1 plus resid 124 and name CA"
pairs["R121_NH2_A124_CA"]="resid 121 and name NH2 plus resid 124 and name CA"
pairs["R121_NH1_V125_CA"]="resid 121 and name NH1 plus resid 125 and name CA"
pairs["R121_NH2_V125_CA"]="resid 121 and name NH2 plus resid 125 and name CA"
pairs["R121_NH1_I213_CA"]="resid 121 and name NH1 plus resid 213 and name CA"
pairs["R121_NH2_I213_CA"]="resid 121 and name NH2 plus resid 213 and name CA"
pairs["R121_NH1_A249_CA"]="resid 121 and name NH1 plus resid 249 and name CA"
pairs["R121_NH2_A249_CA"]="resid 121 and name NH2 plus resid 249 and name CA"


# Loop through each pair and calculate the distance
for key in "${!pairs[@]}"; do
    echo "Calculating distance for $key"
    gmx distance -s conf-md.tpr -f fit.xtc -select "${pairs[$key]}" -oall distance.xvg
    mv distance.xvg "${key}_distance.xvg"
done

echo "All distances calculated and files renamed."
