#!/bin/bash

# Check if the collision_energy is provided as a command line argument
if [ -z "$1" ]; then
    echo "Usage: $0 <collision_energy>"
    exit 1
fi

# Get the collision_energy from the first command line argument
collision_energy="$1"

# List of folders
folders=("AuAu${collision_energy}_00_10" "AuAu${collision_energy}_10_20" "AuAu${collision_energy}_20_30" "AuAu${collision_energy}_30_40" "AuAu${collision_energy}_40_50" "AuAu${collision_energy}_50_60" "AuAu${collision_energy}_60_70" "AuAu${collision_energy}_70_80")

# Loop over each folder
for folder in "${folders[@]}"
do
    scp lpdu318@beluga.computecanada.ca:/lustre03/project/6002853/lpdu318/iEBE-sampler/$folder/event_0/EVENT_RESULTS_MCGlb${folder}_0/spvn_results_MCGlb${folder}_0.h5 /Users/Lipei/Downloads/sampler_test_runs/Beluga_BES_run2
done

