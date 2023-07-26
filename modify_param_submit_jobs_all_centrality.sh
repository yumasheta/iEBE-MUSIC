#!/bin/bash

# Function to process the input file with a given range
process_input_file() {
    python_parameter_file="$1"  # The name of the python file to be modified
    collision_energy="$2"

    for i in {0..7}; do
        modified_centrality_bin=$(printf "%02d_%02d" "$((i*10))" "$((i*10+10))") # Modified centrality bin: 00_10 -> 10_20, etc.

        # Create a temporary directory with a name based on the centrality bin
        iter_temp_dir="temp_${modified_centrality_bin}"
        mkdir -p "$iter_temp_dir"

        # Modify the python_parameter_file using awk and save to a unique temporary file
        awk -v new_pattern="database_name': \"3DMCGlauber_database/MCGlbAuAu${collision_energy}_${modified_centrality_bin}" '{gsub(/database_name.*_[0-9]{2}_[0-9]{2}/, new_pattern)}1' "$python_parameter_file" > "$iter_temp_dir/temp_${modified_centrality_bin}.py"

        ./generate_jobs.py -w "AuAu${collision_energy}_${modified_centrality_bin}" -c McGill -n_urqmd 40 -n_th 40 -par "$iter_temp_dir/temp_${modified_centrality_bin}.py"

        # Remove the temporary directory and its contents
        rm -rf "$iter_temp_dir"

        # Change to the event_0 directory
        cd "AuAu${collision_energy}_${modified_centrality_bin}/event_0"
        
        # Submit the job using sbatch
        sbatch submit_job.pbs
        
        # Go back to the home directory
        cd ../..

    done
}

# Check if the required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <python_parameter_file> <collision_energy>"
    exit 1
fi

python_parameter_file="$1"
collision_energy="$2"

# Call the function with the provided python_parameter_file and collision_energy
process_input_file "$python_parameter_file" "$collision_energy"
