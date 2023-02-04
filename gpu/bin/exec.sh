#!/bin/bash

#SBATCH --partition=t37cluster
#SBATCH --nodelist=phi[28]
#SBATCH --cpus-per-task=1

# Determine number of CPUs to optimize parallel resources
nl=$(lscpu -p | wc -l)
NCPU=$(expr $nl - 4)

echo $NCPU

srun --nodes=1 --ntasks=1 --exclusive --cpus-per-task=$NCPU ./tubenz

