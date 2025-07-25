#!/bin/bash
#SBATCH --job-name=unif_lmm_sims   #Job name	
#SBATCH --mail-type=END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --account=k.ekvall
#SBATCH --qos=k.ekvall-b
#SBATCH --array=1-30
#SBATCH --mail-user=k.ekvall@ufl.edu   # Where to send mail	
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=10:00:00   # Walltime
#SBATCH --output=/blue/k.ekvall/k.ekvall/Simulations/unif_lmm/R1/Output/r_job.%j_%a.out
#Record the time and compute node the job ran on
date; hostname; pwd
#Use modules to load the environment for R
module load R

#Run R script 
Rscript /blue/k.ekvall/k.ekvall/Simulations/unif_lmm/R1/sim_master_R1.R

date