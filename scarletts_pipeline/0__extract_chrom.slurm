#!/bin/bash
#SBATCH --job-name=all
#SBATCH --mail-user=sh700@duke.edu
#SBATCH --mail-type=FAIL
#SBATCH --mem=512M
#SBATCH -c 2
#SBATCH --out=log/extractVcf-1KGP-%j.out 
#SBATCH --error=log/extractVcf-1KGP-%j.out
#SBATCH --time=00:30:00
#SBATCH -p scavenger

set -e
# load modules
module load htslib
# specify tmp path
mkdir -p /scratch/$USER/tmp
mkdir -p /scratch/$USER/cache
export TMPDIR=/scratch/$USER/tmp
export TMP=/scratch/$USER/tmp
export TEMP=/scratch/$USER/tmp
# read inpput
sample=$1
input_vcfgz=$2
output_file=$3
output_header_file=$4

threads=$SLURM_CPUS_ON_NODE

now=$(date +"%T")
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with ${sample} ${input_vcfgz}:"${now}

sample_index=$(
    bgzip -@ $threads -cd $input_vcfgz \
    | awk '$1=="#CHROM" {for(i=1;i<=NF;i++) { if($i=="'$sample'") { print i; } } exit; }'
)

if [ -z $sample_index ]; then
    echo "unable to find sample index for sample $sample"
    exit 1
fi

bgzip -@ $threads -cd $input_vcfgz \
    | grep -v '^#' \
    | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$'$sample_index'}' \
    | awk -v OFS='\t' 'length($4)==1 && length($5)==1 && $7=="PASS"' \
    > $output_file

wc -l $output_file

if [ -n "$output_header_file" ]; then
    echo "generating header file $output_header_file"
    bgzip -@ $threads -cd $input_vcfgz \
        | awk -v OFS='\t' '/^##/ {print $0} /^#/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$'$sample_index'} !/^#/ {exit}' \
        > $output_header_file

    wc -l $output_header_file
fi

now=$(date +"%T")
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> done with ${sample} :"${now}
