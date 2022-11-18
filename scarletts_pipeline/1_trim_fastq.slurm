#!/bin/bash
#SBATCH --job-name=trimFastq
#SBATCH --mail-user=sh700@duke.edu
#SBATCH --mail-type=FAIL
#SBATCH --mem=10G
#SBATCH --time=8:00:00
#SBATCH -c 4
#SBATCH --output=log/step1-trimFastq-%j.out
#SBATCH --error=log/step1-trimFastq-%j.out
#SBATCH -p scavenger

set -e
set -u

############### testing command
# sample="HG00097"
# sbatch ../script/RNAseq-analysis/pipeline_scripts/1_trim_fastq.slurm $sample /hpc/group/allenlab/scarlett/output/RNAseq/1000Genome/$sample/tmp/trimmed_fastq /datacommons/allenlab/scarlett/data/fastq/1000Genome/$sample /hpc/group/allenlab/scarlett/software/Trimmomatic-0.39
###############

# specify tmp path
mkdir -p /scratch/$USER/tmp
mkdir -p /scratch/$USER/cache
export TMPDIR=/scratch/$USER/tmp
export TMP=/scratch/$USER/tmp
export TEMP=/scratch/$USER/tmp
# read input
sample=$1
trimmedfastq_path=$2
rawfastq_path=$3
Trimmomatic_path=$4
real_pileup=$5
# check 
mkdir -p $trimmedfastq_path
tmp_fwd_reads=$trimmedfastq_path/${sample}_FWD_paired.tmp.fq.gz
tmp_rev_reads=$trimmedfastq_path/${sample}_REV_paired.tmp.fq.gz
tmp_fwd_reads_unpaired=$trimmedfastq_path/${sample}_FWD_unpaired.tmp.fq.gz
tmp_rev_reads_unpaired=$trimmedfastq_path/${sample}_REV_unpaired.tmp.fq.gz

fwd_reads=$trimmedfastq_path/${sample}_FWD_paired.fq.gz
rev_reads=$trimmedfastq_path/${sample}_REV_paired.fq.gz
fwd_reads_unpaired=$trimmedfastq_path/${sample}_FWD_unpaired.fq.gz
rev_reads_unpaired=$trimmedfastq_path/${sample}_REV_unpaired.fq.gz

if [ -s "$real_pileup" ]; then
    echo "$real_pileup exists."
    exit
fi

if [ -s "$fwd_reads" \
    -a -s "$rev_reads" \
]; then
    echo "$fwd_reads and $rev_reads exist."
    exit
fi

threads=$SLURM_CPUS_ON_NODE

now=$(date +"%T")
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with ${sample} :"${now}

# use zsh as it waits until process substitions finish before completing.
zsh -c "eval \
java -jar $Trimmomatic_path/trimmomatic-0.39.jar PE -threads $threads -phred33 \
    $rawfastq_path/*1.fastq.gz \
    $rawfastq_path/*2.fastq.gz \
    >(pigz --fast -c > $tmp_fwd_reads) \
    >(pigz --fast -c > $tmp_fwd_reads_unpaired) \
    >(pigz --fast -c > $tmp_rev_reads) \
    >(pigz --fast -c > $tmp_rev_reads_unpaired) \
    ILLUMINACLIP:$Trimmomatic_path/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
"

echo "waiting $(date +%T)"
wait
echo "done waiting $(date +%T)"
# rename completed files

mv $tmp_fwd_reads_unpaired $fwd_reads_unpaired
mv $tmp_rev_reads_unpaired $rev_reads_unpaired
mv $tmp_fwd_reads $fwd_reads
mv $tmp_rev_reads $rev_reads

# echo "fwd_reads lines"
# zcat $fwd_reads | wc -l

now=$(date +"%T")
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> done with ${sample} :"${now}
