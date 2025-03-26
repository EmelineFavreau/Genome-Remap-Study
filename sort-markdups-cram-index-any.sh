#!/bin/bash
##SBATCH -J sortmarkdupcramAny
#SBATCH -A CWALLACE-INTREPID-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --ntasks=38
#SBATCH -t 16:00:00
#SBATCH --mail-type=FAIL
#! %A means slurm job ID and %a means array index
#SBATCH --output=logs/sortmarkdupcramAny_%A_%a.out
#! Errors filename:
#SBATCH --error=logs/sortmarkdupcramAny_%A_%a.err
#!array job, from 1 to 10
#SBATCH --array=1-10

# launch script
# sbatch /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/sort-markdups-cram-index-any.sh BPD-unsorted-bam-list.txt
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment

module load samtools # v1.9
module load picard # 2.9.2

# if a command fails, the shell exit and 
# writes to standard error if a variable is not set and
# writes to standard error a trace of each command
set -eux

# 2 columns: array and sample name
config=$1

# Extract info for the current $SLURM_ARRAY_TASK_ID
# sample=`head -n 1 $config | cut -f 2`
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# 1: not_sorted_bam: eg-bam3-no-sort.bam, the output of the align.sh script
UNS_BAM="/rds/project/rds-hxHuvOdtftk/hpc-work/result/${sample}.GRCh38.bam"

# 2: fasta is available
REF="/rds/project/who1000-1/rds-who1000-cbrc/user/et341/shared/fasta/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# 3: output bam picard-markduplicated, samtools-indexed: eg-bam3.bam
OUT_BAM="/rds/project/rds-hxHuvOdtftk/hpc-work/tmp/${sample}.GRCh38.markduplicates.sorted.bam"

# 4: output metrics from picard: eg-bam3-metrics.txt
OUT_METRICS="/rds/project/rds-hxHuvOdtftk/hpc-work/tmp/${sample}.picard.metrics.txt"

# 5: output counts from samtools: eg-bam3-counts.txt
OUT_COUNTS="/rds/project/rds-hxHuvOdtftk/hpc-work/tmp/${sample}.markduplicates.sorted.counts.before.after.txt"

# 6: temporary_sorted_bam from picard: eg-bam3-sort.bam
SORTED_TMP="/rds/project/rds-hxHuvOdtftk/hpc-work/tmp/${sample}.picard.sorted.tmp.bam"

# 7: lots of space tempory folder for picard: tmp
TMP_DIR="/rds/project/rds-hxHuvOdtftk/hpc-work/tmp"

# 8: output cram picard-markduplicated, samtools-indexed, no compression: eg-bam3.cram
OUT_CRAM="/rds/project/rds-yWykMIGjpDE/WgsReMap2024/alignments/${sample}.GRCh38.markduplicates.sorted.cram"

# count reads before sort/markduplicate
samtools view -c ${UNS_BAM} > ${OUT_COUNTS}

# sort reads 
samtools sort -o ${SORTED_TMP} ${UNS_BAM}

# remove the unsorted bam
rm -rf ${UNS_BAM}

# takes unsorted bam, markduplicates, output as bam
java -jar $PICARD_HOME/picard.jar MarkDuplicates I=${SORTED_TMP}  ASSUME_SORT_ORDER=coordinate M=${OUT_METRICS} R=${REF} O=${OUT_BAM} TMP_DIR=${TMP_DIR} 

# take a bam, output a cram
samtools view -T ${REF} -C -o ${OUT_CRAM} ${OUT_BAM}

# create an index file
samtools index ${OUT_CRAM}

# remove the temporary sorted bam and bams
rm -rf ${SORTED_TMP}
rm -rf ${OUT_BAM}

# counts reads after 
samtools view -T ${REF} -c ${OUT_CRAM} >> ${OUT_COUNTS}
