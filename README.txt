## From GRCh37 bam to GRCh38 cram

1. **Send mapping jobs** 
`bash send-align-sort-jobs-any.sh ${fileLines} ${alignConfig} ${sortConfig}`
*e.g.* `bash send-align-sort-jobs-any.sh 1-10 array-sample-path.txt array-sample.txt`
This command takes three arguments: (1) the array for Slurm, (2) a text file with bam paths (3 columns: array, sample name, bam path), (3) a text file with sample names (2 columns: array, sample name). It will send in the queue two jobs:
	1. for the script `align-any.sh`: array job slurm script to run `bazam`, `bwa`, `samtools` to create a temporary GRCh38 bam file;
    2. for the script `sort-markdups-cram-index-any.sh`:   array job slurm script to run `samtools` and `picard` to create final sorted and mark-duplicated CRAM files.                  

2. **Send variant calling jobs**
`bash send-deepvariant-jobs.sh ${fileLines}`
*e.g.* `bash send-deepvariant-jobs.sh 1-10`

This command takes one argument (the array for Slurm). It will send in the queue a job for the script `deepvariant-calls-variants.sh`. This is an array job slurm script to run nextflow nf-core sarek pipeline to create annotated VCF files using Deepvariant, VEP, SnpEff.
