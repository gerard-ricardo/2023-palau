#!/bin/bash
#SBATCH --job-name=myjob
#SBATCH --account=a_msel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=2:00:00
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=g.ricardo@uq.edu.au


# Generate a timestamp
#timestamp=$(date +%Y%m%d-%H%M%S)

#cp "./R_scripts/1_Bunya_master.R" "./script_txts/1_Bunya_master1.R"


module load r/4.2.1-foss-2022a



srun Rscript -e install.packages('dartR', repos="https://cran.csiro.au/")