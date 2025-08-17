#!/bin/bash
#SBATCH --job-name="Electro"
#SBATCH --nodelist=gpu-nvidia
#SBATCH --mem=120G # Memory needed per task
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks=32                # Maximum number of MPI tasks
#SBATCH --time=12:00:00            # Job time limit (hh:mm:ss)


# load mpi
module load openmpi
export OMPI_MCA_psec=^munge

echo "N_THREADS,REP1,REP2,REP3"
for i in 1 2 4 8 16 32; do
    echo -n "$i"
    for r in $(seq 1 3); do
        export OMP_NUM_THREADS=$i
        result=$(./electromain TraxStick3 0)
        echo -n ",$result"
    done
    echo
done


echo "All runs completed!"
