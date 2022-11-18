#!/bin/bash
#SBATCH --job-name=RNAseq_alignment
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --mail-user=sh700@duke.edu
#SBATCH --mail-type=FAIL
#SBATCH --out=log/step2-RNAseqAlignment-%j.out
#SBATCH --error=log/step2-RNAseqAlignment-%j.out
#SBATCH -p scavenger
#SBATCH --time=12:00:00

############### testing command
# sbatch ../script/RNAseq-analysis/pipeline_scripts/2__star_alignment.slurm HG00097 /datacommons/allenlab/scarlett/reference/hg19 /datacommons/allenlab/scarlett/reference/STARIndex /datacommons/allenlab/hg19/filter/gencode.v19.annotation.level12.gtf /hpc/group/allenlab/scarlett/output/RNAseq/1000Genome/HG00097/tmp/trimmed_fastq /hpc/group/allenlab/scarlett/output/RNAseq/1000Genome/HG00097/tmp/star_2pass_WASP_m10 /hpc/group/allenlab/scarlett/output/RNAseq/1000Genome/HG00097/HG00097.no_chr.content.SNPs.hets.vcf.gz 10
###############

# cmd to generate star index
#srun STAR --runMode genomeGenerate --genomeDir $star_ind --genomeFastaFiles $ref/hg19.fa --sjdbGTFfile $Annotation --runThreadN 8 --outFileNamePrefix $star_ind


set -e
set -u

# load modules
module load samtools
module load STAR/2.7.5c
# specify tmp path
mkdir -p /scratch/$USER/tmp
mkdir -p /scratch/$USER/cache
export TMPDIR=/scratch/$USER/tmp
export TMP=/scratch/$USER/tmp
export TEMP=/scratch/$USER/tmp
# read input
sample=$1
ref=$2
star_ind=$3
Annotation=$4
fastq_path=$5
alignedbam_path=$6
vcfgz=$7
N=$8
real_pileup=$9

threads=$SLURM_CPUS_ON_NODE

in_fwd_fastq=$fastq_path/${sample}_FWD_paired.fq.gz
in_rev_fastq=$fastq_path/${sample}_REV_paired.fq.gz
out_bam=$alignedbam_path/Aligned.sortedByCoord.out.bam

if [ "$out_bam" -nt "$in_fwd_fastq" -a -s "$out_bam" ]; then
    echo "$out_bam exists and non-empty, newer than $in_fwd_fastq."
    exit
fi

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with ${sample} $(date +%T)"
mkdir -p $alignedbam_path

STAR --twopassMode Basic \
--runThreadN $threads \
--genomeDir $star_ind \
--readFilesIn $in_fwd_fastq $in_rev_fastq \
--alignEndsType EndToEnd \
--waspOutputMode SAMtag \
--varVCFfile <(zcat $vcfgz | grep -E "^[^#]" | sed -e 's/^/chr/') \
--outSAMattributes vA vG vW \
--outFilterMismatchNmax $N \
--outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx \
--readFilesCommand "zcat" \
--limitBAMsortRAM 107374182400 \
--outFileNamePrefix $alignedbam_path/
echo "Finish STAR 2pass EndtoEnd alignment"

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> indexing with ${sample} $(date +%T)"

srun samtools index -@ $threads $out_bam

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> done with ${sample} $(date +%T)"
