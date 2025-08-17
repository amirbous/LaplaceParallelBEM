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

echo "N_MPI,N_OMP,REP1s,REP2s,REP3s,REP4s,REP5s"

for i in 1 2 4 8 16 32; do
    max_j=$((64 / i))   # maximum threads allowed for this i
    j=1
    while [ $j -le $max_j ]; do
        echo -n "$i,$j"
        for r in $(seq 1 5); do
            export OMP_NUM_THREADS=$j
            result=$(srun -n $i ./electroparal Cylinder_5)
            echo -n ",$result"
        done
        echo
        j=$((j * 2))  # power of two progression
    done
done



echo "All runs completed!"
