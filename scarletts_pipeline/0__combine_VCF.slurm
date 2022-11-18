#!/bin/bash
#SBATCH --get-user-env
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH --mail-user=sh700@duke.edu
#SBATCH --mail-type=FAIL
#SBATCH --out=combine-vcf-%j.out
#SBATCH --error=combine-vcf-%j.out
#SBATCH -p scavenger
#SBATCH --time=00:30:00

set -u
set -e

# load module
module load htslib
# specify tmp path
mkdir -p /scratch/$USER/tmp
mkdir -p /scratch/$USER/cache
export TMPDIR=/scratch/$USER/tmp
export TMP=/scratch/$USER/tmp
export TEMP=/scratch/$USER/tmp

# read input
sample=$1
tmp_dir=$2
bi_vcfgz=$3
bihet_vcfgz=$4

now=$(date +"%T")
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with ${sample} :"${now}

cat $tmp_dir/${sample}_chr{1..22}.vcf > $tmp_dir/${sample}_allchr.vcf
cat $tmp_dir/${sample}_allchr.vcf | awk '{ if ($10=="1|0" || $10=="0|1") {print} }' > $tmp_dir/${sample}_allchr.hets.vcf

######################################### all SNPs
cat $tmp_dir/${sample}_header.vcf \
    $tmp_dir/${sample}_allchr.vcf \
    > $tmp_dir/${sample}.no_chr.content.SNPs.vcf

bgzip -@2 -c -f $tmp_dir/${sample}.no_chr.content.SNPs.vcf > $bi_vcfgz \
    && tabix -fp vcf $bi_vcfgz

if [ -s "$bi_vcfgz" ]; then
    rm $tmp_dir/${sample}_allchr.vcf
    rm $tmp_dir/${sample}.no_chr.content.SNPs.vcf
fi

######################################### hets
cat $tmp_dir/${sample}_header.vcf \
    $tmp_dir/${sample}_allchr.hets.vcf \
    > $tmp_dir/${sample}.no_chr.content.SNPs.hets.vcf

bgzip -@2 -c -f $tmp_dir/${sample}.no_chr.content.SNPs.hets.vcf > $bihet_vcfgz \
    && tabix -fp vcf $bihet_vcfgz

if [ -s "$bihet_vcfgz" ]; then
    rm $tmp_dir/${sample}_allchr.hets.vcf
    rm $tmp_dir/${sample}.no_chr.content.SNPs.hets.vcf
    rm $tmp_dir/${sample}_chr*.vcf
fi

now=$(date +"%T")
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> done with ${sample} :"${now}
