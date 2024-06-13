#!/bin/bash
#SBATCH --job-name=fastq
#SBATCH --mail-type=ALL # Valid values: BEGIN, END, FAIL, REQUEUE and ALL.
#SBATCH --mail-user=baponterolon@tulane.edu
#SBATCH --output=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC/fastQC_files/fastqc.out
#SBATCH --error=/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC/fastQC_files/fastqc.err
#SBATCH --time=0-24:00:00
#SBATCH --partition=centos7
#SBATCH --cpus-per-task=20  #: Cpus per Task
#SBATCH --nodes=1            #: Number of Nodes


#### load modules###
module load fastqc
module load gnuparallel
module load anaconda3/2023.07
#unset PYTHONPATH


conda activate py310 # This is to load a virtual environment (py310) from which to load 'multiqc'
###########################################################################################
<<Note_toUser

Modified by BAR 2024-06-12
The script below serves to perform most of the QC steps needed before aligning and 
variant calling.
In a previous version the SLURM header ran as an array job. This is not necessary because we are running a `for` loop in the script.

Softlinks to files only need to be made once and it's best not to work with the original 
fastq or fastq.gz files. 

FastQC will generate QC reports for each sequence file. This step is very slow to run.
MultiQC will then take the fastQC output files and compile them into a single QC report.
Note_toUser
###########################################################################################

echo "Start Job"

### Create folders fastQC_files and multiqc in 2_fastQC FIRST
### set working directory (file output directory) ###
### Variables for the fastq working directory and the input files directory ###
WD="/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC/"
INPUT="/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/1_data/"
cd $WD

### creating softlinks for fastq files to current directory ###
ln -s =/lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/1_data/* .

### running FastQC ###
echo "Start FastQC"
for fastq in /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC/*;do
fastqc $fastq -o /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC/fastQC_files/
done
echo "End FastQC"

### running MultiQC (using Fastqc files) ###
echo "Start MultiQC"
multiqc . -o /lustre/project/svanbael/bolivar/Mimulus_sequences/mim3_bioinformatics/ddRAD/2_fastQC/multiqc
echo "End MultiQC"

### Unsetting environmentc ###
conda deactivate
module purge #unloads all module

echo "End Job"
