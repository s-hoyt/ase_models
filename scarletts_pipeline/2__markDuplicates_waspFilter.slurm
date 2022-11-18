#!/bin/bash
#SBATCH --job-name=markDuplicates
#SBATCH --mem=64G
#SBATCH -c 6
#SBATCH --mail-user=sh700@duke.edu
#SBATCH --mail-type=FAIL
#SBATCH --out=log/step2-markduplicates-WASP-%j.out
#SBATCH --error=log/step2-markduplicates-WASP-%j.out
#SBATCH -p scavenger
#SBATCH --time=8:00:00

############### testing command
# sbatch ../script/RNAseq-analysis/pipeline_scripts/2__markDuplicates_waspFilter.slurm /hpc/group/allenlab/scarlett/output/RNAseq/1000Genome/HG00097/tmp/star_2pass_WASP_m10
###############

set -e
set -u

hostname 2>&1
# load modules
module load samtools
module load Picard/2.18.2
# specify tmp path
mkdir -p /scratch/$USER/tmp
mkdir -p /scratch/$USER/cache
export TMPDIR=/scratch/$USER/tmp
export TMP=/scratch/$USER/tmp
export TEMP=/scratch/$USER/tmp
# read input
alignedbam_path=$1
real_pileup=$2

in_bam=$alignedbam_path/Aligned.sortedByCoord.out.bam
out_bam=$alignedbam_path/Aligned.sortedByCoord.out.picard_markdup.WASP_filtered.bam

if [ -s "$real_pileup" ]; then
    echo "$real_pileup exists."
    exit
fi
if [ "$out_bam" -nt "$in_bam" \
    -a -s "$out_bam" \
    -a "${out_bam}.bai" -nt "$out_bam" \
]; then
    echo "$out_bam exists and non-empty, newer than $in_bam."
    exit
fi

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start $(date +%T)"

java -Xmx48G \
    -jar $PICARD MarkDuplicates \
    I=$in_bam \
    O=$alignedbam_path/Aligned.sortedByCoord.out.picard_markdup.bam \
    M=$alignedbam_path/picard.marked_dup_metrics.txt

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> indexing $(date +%T)"

# srun samtools index -@ 4 $alignedbam_path/Aligned.sortedByCoord.out.picard_markdup.bam

samtools view -@ 4 -h $alignedbam_path/Aligned.sortedByCoord.out.picard_markdup.bam \
    -e '!([vW] >= 2)' -b -o $out_bam
srun samtools index -@ 4 $out_bam

# echo "Finish filtering"
# echo "Start samtools flagstat for picard"
# srun samtools flagstat $alignedbam_path/Aligned.sortedByCoord.out.picard_markdup.WASP_filtered.bam > $alignedbam_path/flagstat_markdup_picard.txt
# echo "Start counting unique read pair"
# srun samtools view -@ 4 -bh -q 30 -f 3 -F 2316 -c $alignedbam_path/Aligned.sortedByCoord.out.picard_markdup.WASP_filtered.bam
# srun samtools view -@ 4 -bh -q 30 -f 3 -F 3340 -c $alignedbam_path/Aligned.sortedByCoord.out.picard_markdup.WASP_filtered.bam

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> done $(date +%T)"
