#!/bin/bash

# This master script sends a slurm job to obtain deepvariant vcf

# launch script to obtain, for instance 238 (those files listed on lines 11 to 249)
# bash /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/send-deepvariant-jobs.sh 11-249
# $1 is the line number of the first sample, to have a VCF file created, and the last one

fileLines=${1}

### task: call variants with deepvariant
# fileLines=`cat this-vcf-list | tr '\n' ',' | sed 's/,$/%14\n/'`
# update the file lines
sed -i "s/1-10/${fileLines}/g" /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/deepvariant-calls-variants.sh

# send the slurm array script 
sbatch /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/deepvariant-calls-variants.sh

### tidy space
sed -i "s/${fileLines}/1-10/g" /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/deepvariant-calls-variants.sh
