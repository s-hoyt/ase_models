#!/bin/bash
#SBATCH --job-name=submit-pipelines-1000Genome-GBR
#SBATCH --mem=512M
#SBATCH -c 1
#SBATCH --mail-type=END,FAIL
#SBATCH -p scavenger
#SBATCH --oversubscribe
#SBATCH --no-requeue
#SBATCH --dependency=singleton
#SBATCH --out=log/submit-pipelines-%j.out
#SBATCH --error=log/submit-pipelines-%j.out

module load Parallel

set -x
set -u
set -e
########### testing commands
# sbatch ../script/RNAseq-analysis/pipeline_scripts/submit_pipelines.slurm
###########

###########  1000 Genome samples
path=/datacommons/allenlab/scarlett/data/fastq/1000Genome
#file=$path/British.txt GBR #188/2=94 #8-86-86 HG00104 HG00124 HG00134 HG00135 HG00152 HG00156 HG00247 HG00249
#file=$path/Finnish.txt FIN #190/2=95 #3-92-92 HG00312 HG00359 HG00377
#file=$path/Tuscan.txt  TSI #186/2=93 #2-91-91 NA20537 NA20816
#file=$path/Utah.txt    CEU #182/2=91 #2-89-89 NA11993 NA07346
#file=$path/Yoruba.txt  YRI #178/2=89 #2-87-87 NA18487 NA19150
# 94+95+93+91+89 = 465 - (9+3+3+2+2) = 445 | 445

ancestry="British"
file=$path/${ancestry}.txt
start=1
count=100
samples=$(cat ${file} | awk '{print $1}' | uniq | tail -n+${start} | head -n $count)

function submit_and_wait_sample() {
    sample=$1
    ancestry="GBR"
    #### 1000 Genome
    base_working_dir=/hpc/group/allenlab/stephanie/output/RNAseq/1000Genome
    vcf_dir=/datacommons/allenlab/scarlett/data/VCF/1000_genome/20130502/bgzip
    raw_fq_dir=/datacommons/allenlab/scarlett/data/fastq/1000Genome
    read_length=75
    simulation_depth=100
    #
    job_id=$(sbatch \
        --parsable \
        --job-name=${sample}-${ancestry}-pipeline \
        /hpc/home/sh700/ase_models/scarletts_pipeline/0_extract_vcf.slurm \
        --sample $sample \
        --vcf-dir $vcf_dir \
        --raw-fq-dir $raw_fq_dir \
#        --read-length $read_length \
#        --simulation-depth $simulation_depth \
#        --base-working-dir $base_working_dir \
#        --is-1000Genome \
#        --simulate-allSNPs \
	--force
    )
    ../script/RNAseq-analysis/pipeline_scripts/wait_for_jobs.sh $job_id
}
export -f submit_and_wait_sample
# options for pipeline
# --is-1000Genome    --> default: individual sample
# --simulate-allSNPs --> default: only mutate het SNPs 
# --keep-tmp         --> default: remove tmp directory
# --random           --> default: even haplotype mode so that p = 0.5
# --force            --> default: check success #file
#samples="HG00130 HG00132" #HG00265 HG00338 HG00341 HG00342 HG00343 HG00344 HG00345 HG00346 HG00349 HG00353 HG00355 HG00356 HG00361 HG00362 HG00365 HG00366 HG00367 HG00369 HG00371 HG00372 HG00373 HG00375 HG00376 HG00378 HG00379 HG00380 HG00381 HG00382 HG00383 HG00384 HG01790 NA12873 NA12874 NA12889 NA12890 NA18488 NA18868 NA19117 NA20581 NA20585 NA20758 NA20783 NA20785 NA20786 NA20796 NA20797 NA20799 NA20800 NA20803 NA20807 NA20809 NA20810 NA20811 NA20812 NA20813 NA20814 NA20815 NA20819 NA20826 NA20828"
parallel \
    -j 5 \
    --delay 1 \
    -v --line-buffer \
    submit_and_wait_sample \
    ::: $samples

#start 20 samples,submit each one after 1s of submitting the previous one
