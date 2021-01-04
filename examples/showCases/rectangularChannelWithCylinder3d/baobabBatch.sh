#!/bin/sh

#SBATCH --partition=cui-EL7,parallel-EL7,mono-EL7

#SBATCH --time=4-00:00:00

#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1

#SBATCH --job-name=kamil
#SBATCH --output=slurm-%J.out

# Have the same in the .bashrc
module load CMake
module load foss/2019a
module load ImageMagick

srun ./rectangularChannelWithCylinder3d 70 20 6 200 200 200 1 1 0.6