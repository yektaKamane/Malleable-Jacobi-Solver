#!/bin/bash

# Check if required SLURM environment variables exist
if [[ -z "$SLURM_NNODES" || -z "$SLURM_NODELIST" || -z "$SLURM_NTASKS_PER_NODE" || -z "$SLURM_NPROCS" ]]; then
    echo "Error: One or more SLURM environment variables are missing."
    exit 1
fi

# Expand the node list (e.g., "awnode[05-06]" -> "awnode05 awnode06")
NODELIST_EXPANDED=$(scontrol show hostnames $SLURM_NODELIST)

# Create or overwrite the nodes file
> mynodelist

# Loop through each node and write it to the file SLURM_NTASKS_PER_NODE times
for node in $NODELIST_EXPANDED; do
    for ((i = 0; i < SLURM_NTASKS_PER_NODE; i++)); do
        echo "host $node" >> mynodelist
    done
done

echo "nodes file created successfully!"

# Calculate the total number of non-empty lines in mynodelist
# new_node_number=$(grep -cve '^\s*$' mynodelist)

# If exit code is 0, exit
./charmrun +p$SLURM_NPROCS ./jacobi 1 1 $SLURM_NPROCS 100 +balancer GreedyLB ++nodelist ./mynodelist ++server ++server-port 1234 ++verbose

# If exit code is 102
while [[ $? -eq 102 ]]; do
    
    new_node_number=$(grep -c '^host' mynodelist)
    echo "Jacobi exited with code 102. Restarting with $new_node_number nodes."

    ./charmrun +p$new_node_number ./jacobi 1 1 $new_node_number 100 +balancer GreedyLB ++nodelist ./mynodelist ++server ++server-port 1234 ++verbose +restart log
    exit_code=$?
    if [[ $exit_code -eq 0 ]]; then
        echo "Jacobi completed successfully after restart. Exiting."
        exit 0
    fi
done
