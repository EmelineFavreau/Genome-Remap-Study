#!/bin/bash

#SBATCH -J alignAny
#SBATCH -A CWALLACE-INTREPID-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=16:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=logs/alignAny_%A_%a.out
#! %A means slurm job ID and %a means array index
#! Errors filename:
#SBATCH --error=logs/alignAny_%A_%a.err
#!array job, from 1 to 10
#SBATCH --array=1-10


# launch script
# sbatch /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/align-any.sh BPD-bam-list.txt this-list-BPD BPD-bam-list.txt
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment

module load samtools # v1.9
module load bwa-0.7.17-gcc-5.4.0-42mry2g

# if a command fails, the shell exit and 
# writes to standard error if a variable is not set and
# writes to standard error a trace of each command
set -eux

# list of bam paths and sample name
config=$1

# Extract info for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

input_bam=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

out_bam="/rds/project/rds-hxHuvOdtftk/hpc-work/result/${sample}.GRCh38.bam"

# count reads before re-aligning
samtools view -c ${input_bam} > /rds/project/rds-hxHuvOdtftk/hpc-work/tmp/${sample}-bam-counts.txt

# remove alignment information, change from bam to fastq, align to GRCh38, save output as bam
java -jar /rds/user/ef479/hpc-work/bazam.jar -bam ${input_bam} | bwa mem -Y -K 100000000 -t 15 -p /rds/project/who1000-1/rds-who1000-cbrc/user/et341/shared/fasta/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa - | samtools view -bS -o ${out_bam}

# count reads after re-aligning
#echo "did not count" >> /rds/project/rds-hxHuvOdtftk/hpc-work/tmp/${sample}-bam-counts.txt
samtools view -c ${out_bam} >> /rds/project/rds-hxHuvOdtftk/hpc-work/tmp/${sample}-bam-counts.txt
