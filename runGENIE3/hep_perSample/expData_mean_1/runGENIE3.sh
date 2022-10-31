#!/bin/bash

#SBATCH --job-name=GENIE3
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --time=480:00:00
#SBATCH --mail-user=tugba.agaoglu@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/tagaoglu/output/output_SigMat_%j.o
#SBATCH --error=/data/users/tagaoglu/error/error_SigMat_%j.e
#SBATCH --array=0-1

WORK_DIR=/data/projects/p728_scRNA_HCC/RSTUDIO/analysis/tagaoglu/runGENIE3/hep_perSample/expData_mean_1

cd $WORK_DIR

SCRIPTS=(*.R)

#To run R commands without starting an R session
Rscript ${SCRIPTS[$SLURM_ARRAY_TASK_ID]}
