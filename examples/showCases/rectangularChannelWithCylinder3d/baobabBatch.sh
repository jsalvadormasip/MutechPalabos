#!/bin/sh

#SBATCH --partition=cui-EL7,parallel-EL7,mono-EL7,shared-EL7,mono-shared-EL7

#SBATCH --time=0-12:00:00

#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1

#SBATCH --job-name=kamil
#SBATCH --output=slurm-%J.out

# SBATCH --array=20,50,100,150
# srun ./rectangularChannel3dWithCylinder3d 100 40 6 ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_TASK_ID} 150 1 1 0.6

# Have the same in the .bashrc
module load CMake
module load foss/2019a
module load ImageMagick

srun ./rectangularChannel3dWithCylinder3d 1200 40 6 21 6 75 1 1 0.6