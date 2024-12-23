#!/bin/bash
#SBATCH -p xhacnormalb
#SBATCH -N 1
#SBATCH --ntasks-per-node=80
#SBATCH --cpus-per-task=1
##SBATCH --exclusive

export PATH=/work/share/acumivbyb3/matlab2023a/bin:$PATH

matlab -nodesktop -nosplash -nodisplay -r NESS_FOUR_n1nL_EQ