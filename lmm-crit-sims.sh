#!/bin/bash
#SBATCH --job-name=lmm-crit-sims   #Job name	
#SBATCH --mail-type=END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --account=
#SBATCH --qos=
#SBATCH --array=1-30
#SBATCH --mail-user=   # Where to send mail	
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=10:00:00   # Walltime
#SBATCH --output=
#Record the time and compute node the job ran on
date; hostname; pwd
#Use modules to load the environment for R
module load R

#Run R script 
Rscript simulation_master.R

date