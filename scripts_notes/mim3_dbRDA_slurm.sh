#!/bin/bash

#SBATCH --qos=normal
#SBATCH --job-name MIM3_dbRDA_model
#SBATCH --error /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/MIM3_dbRDA_models.error     
#SBATCH --output /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/MIM3_dbRDA_models.output  
#SBATCH --time=23:00:00
#SBATCH --mem=128000 #Up to 256000, the maximum increases queue time
#SBATCH --nodes=2               # One is enough. When running MPIs, anywhere from 2-4 should be good.
#SBATCH --ntasks-per-node=1    # Number of Tasks per Node (MPI processes)
#SBATCH --cpus-per-task=20      # Number of cpu-cores (threads) per task (OMP threads)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=baponterolon@tulane.edu
#SBATCH --partition=centos7    #This is important to run the latest software versions

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK #$SLURM_CPUS_PER_TASK stores whatever value we assign to cpus-per-task, and is therefore our candidate for passing to OMP_NUM_THREADS

##Modules/Singularity
module load anaconda3/2020.07
source activate virtual_env

# Basic session info
echo Start Job

echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

#Run R script
R CMD BATCH mim3_dbRDA_models.R # Add your script here

module purge

echo End Job
