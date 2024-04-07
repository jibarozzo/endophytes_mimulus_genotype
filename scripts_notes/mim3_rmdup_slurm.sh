#!/bin/bash

#SBATCH --qos=long
#SBATCH --job-name MIM3_ddRAD_rmdup
#SBATCH --error /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/MIM3_ddRAD_rmdup.error     
#SBATCH --output /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/MIM3_ddRAD_rmdup.output  
#SBATCH --time=2-00:00:00
#SBATCH --mem=128000 #Up to 256000, the maximum increases queue time
#SBATCH --nodes=1               # One is enough. When running MPIs, anywhere from 2-4 should be good.
#SBATCH --ntasks-per-node=1    # Number of Tasks per Node (MPI processes)
#SBATCH --cpus-per-task=20      # Number of cpu-cores (threads) per task (OMP threads)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=baponterolon@tulane.edu

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK # Set OMP_NUM_THREADS to the number of CPUs per task we asked for.

##Modules/Singularity
module load anaconda

# Basic session info
echo Start Job

echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

#Run Python script
filePrefix=("BAR24" "BAR26" "BAR28" "BAR29" "BAR31" "BAR32" "BAR33" "BAR34" "BAR35" "BAR36" "BAR42")

for prefix in "${filePrefix[@]}"; do
    python rmdup_molbarcodes.py -p "$prefix" -s fastq.gz
done
module purge

echo End Job
