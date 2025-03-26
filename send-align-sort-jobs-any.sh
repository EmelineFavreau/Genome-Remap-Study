#!/bin/bash

# This master script sends a slurm job to align bam files to v38
# then as soon as the array job is assigned a JobID
# the script creates a second slurm job script to sort the produced bam

# launch script to obtain, for instance 238 sorted mapped bam files (those files listed on lines 11 to 249)
# bash /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/send-align-sort-jobs-any.sh 11-249 BPD-bam-list.txt BPD-unsorted-bam-list.txt
# $1 is the line number of the first bam file to realign and the last one
# $2 is the config file for the first script
# $3 is the config file for the second script

fileLines=${1}

### task 1: align bam files to v38
# Specify the path to the config file: first column: arrayID, second column name of sample, third column: full path to input bam
# update the file lines
sed -i "s/1-10/${fileLines}/g" /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/align-any.sh
# send the slurm array script and capture the JOB ID into a variable
jobID=$(sbatch --parsable /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/align-any.sh ${2})



### task 2: sort and markduplicate the v38 bam 
# Specify the path to the config file: first column is arrayID, second column is name of sample
# update the file lines
sed -i "s/1-10/${fileLines}/g" /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/sort-markdups-cram-index-any.sh
# send the slurm array script once all jobs are terminated successfully
jobID2=$(sbatch --dependency=aftercorr:${jobID} --parsable /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/sort-markdups-cram-index-any.sh ${3})


### task 3: monitor efficacity
sbatch --dependency=aftercorr:${jobID2} /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/append-log-with-seff.sh $jobID2 sortmarkdupcramAny

### tidy space
sed -i "s/${fileLines}/1-10/g" /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/align-any.sh
sed -i "s/${fileLines}/1-10/g" /rds/project/rds-hxHuvOdtftk/hpc-work/scripts/sort-markdups-cram-index-any.sh
