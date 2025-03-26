#!/bin/bash
#SBATCH -J deepvariantCallVariants
#SBATCH -A CWALLACE-INTREPID-SL2-CPU
#SBATCH --partition icelake-himem
#SBATCH --nodes=1
#SBATCH --time=12:00:00 
#SBATCH --ntasks=38
#SBATCH --mail-type=FAIL
#! %A means slurm job ID
#SBATCH --output=logs/deepvariantCallVariants_%A_%a.out
#! Errors filename:
#SBATCH --error=logs/deepvariantCallVariants_%A_%a.err
#!array job, from 1 to 10
#SBATCH --array=1-10

# launch script: sbatch /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/deepvariant-calls-variants.sh
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment
module load singularity

# parameters for Nextflow
export JAVA_HOME=/home/ef479/bin/jdk-20.0.1
export NXF_SINGULARITY_CACHEDIR=/rds/user/ef479/hpc-work/nextflow-singularity-cache
export TOWER_ACCESS_TOKEN=eyJ0aWQiOiA4NjUxfS5lNDc1NDZmOGI2Y2FmMGY3N2NiNzc4ODI0NDEwZWEzMGI2NTQ2ZGRi
export NXF_TEMP=/rds/project/rds-yWykMIGjpDE/WgsReMap2024/tmp
export NXF_OPTS='-Xms1g -Xmx4g'

# list of bam paths and sample name
config="/rds/project/rds-hxHuvOdtftk/hpc-work/PID-bam-list.txt"

# Extract info for the current $SLURM_ARRAY_TASK_ID (e.g. W000288_A)
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# lane: either the character after _; or if not present; assign A
# Lane ID, used when the sample is multiplexed on several lanes. 
lane=$(echo ${sample} | sed "s/^.*_//g")
[ -z "$lane" ] && lane="A"

# cram crai
input_cram="/rds/project/rds-yWykMIGjpDE/WgsReMap2024/alignments/${sample}.GRCh38.markduplicates.sorted.cram"
input_crai="/rds/project/rds-yWykMIGjpDE/WgsReMap2024/alignments/${sample}.GRCh38.markduplicates.sorted.cram.crai"

# make samplesheet
rm -rf /rds/project/rds-hxHuvOdtftk/hpc-work/sarek-test/input/${sample}.samplesheet.csv
echo "patient,sample,lane,cram,crai" > /rds/project/rds-hxHuvOdtftk/hpc-work/sarek-test/input/${sample}.samplesheet.csv
echo "${sample},${sample},${lane},${input_cram},${input_crai}" >> /rds/project/rds-hxHuvOdtftk/hpc-work/sarek-test/input/${sample}.samplesheet.csv


# Intervals are parts of the chopped up genome used to speed up preprocessing
# default is GATK v38
# --nucleotides_per_second 2000000 parallelize variant calling across genomic chunks of roughly similar sizes
rm -rf /rds/project/rds-yWykMIGjpDE/WgsReMap2024/VCF/${sample}-deepvariant-v38
nextflow -c /rds/project/rds-hxHuvOdtftk/hpc-work/slurm_for_nf.conf run nf-core/sarek --input sarek-test/input/${sample}.samplesheet.csv --outdir /rds/project/rds-yWykMIGjpDE/WgsReMap2024/VCF/${sample}-deepvariant-v38 -with-tower --tools deepvariant,snpeff,vep,merge --step variant_calling --nucleotides_per_second 2000000 -r 3.4.1 --cleanup true

