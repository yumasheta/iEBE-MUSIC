#!/bin/bash

# Check if the collision_energy is provided as a command line argument
if [ -z "$1" ]; then
    echo "Usage: $0 <collision_energy>"
    exit 1
fi

# Get the collision_energy from the first command line argument
collision_energy="$1"

# List of folders
folders=("AuAu${collision_energy}_0_10" "AuAu${collision_energy}_10_20" "AuAu${collision_energy}_20_30" "AuAu${collision_energy}_30_40" "AuAu${collision_energy}_40_50" "AuAu${collision_energy}_50_60" "AuAu${collision_energy}_60_70" "AuAu${collision_energy}_70_80")

# Loop over each folder
for folder in "${folders[@]}"
do
    # Change to the event_0 directory
    cd "$folder/event_0"
    
    # Submit the job using sbatch
    sbatch submit_job.pbs
    
    # Go back to the home directory
    cd ../..
done

