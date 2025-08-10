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

export OMP_NUM_THREADS=64


for i in $(seq 1 32); do
    for r in $(seq 1 3); do
        echo "######################################"
        echo "running with $i cores!, Rep $r"
        echo "######################################"
        srun -n $i ./electroparal Cylinder_8
    done
done
echo "All runs completed!"

