#!/bin/bash
#SBATCH --get-user-env
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --mail-user=sh700@duke.edu
#SBATCH --mail-type=FAIL
#SBATCH --out=log/extract-vcf-%j.out
#SBATCH --error=log/extract-vcf-%j.out
#SBATCH -p scavenger
#SBATCH --time=02:00:00
#SBATCH --oversubscribe
#SBATCH --no-requeue

############### testing command
# sbatch -p common ../script/RNAseq-analysis/pipeline_scripts/0_extract_vcf.slurm NA06989 /hpc/group/allenlab/scarlett/output/RNAseq/NA06989/tmp/VCF /datacommons/allenlab/scarlett/data/VCF/1000_genome/20130502/bgzip/
###############

set -e
set -u

function get_script_dir() {
    local script_path
    if [ -n "$SLURM_JOB_ID" ] ; then
        script_path=$(realpath $(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}' | cut -d' ' -f1 ))
    else
        script_path=$(realpath $0)
    fi
    echo $(dirname "$script_path")
}
scripts_dir=$( get_script_dir )

# specify tmp path
mkdir -p /scratch/$USER/tmp
mkdir -p /scratch/$USER/cache
export TMPDIR=/scratch/$USER/tmp
export TMP=/scratch/$USER/tmp
export TEMP=/scratch/$USER/tmp
# read input
sample=$1
tmp_dir=$2
vcf_dir=$3
bi_vcfgz=$4
bihet_vcfgz=$5
filtered_bi_vcfgz=$6
filtered_bihet_vcfgz=$7
genotypingError_sites=$8
real_pileup=$9

if [ -s "$filtered_bihet_vcfgz" -a -s "$filtered_bi_vcfgz" ]; then
    if [ $(zcat "$filtered_bihet_vcfgz" | wc -l) -gt 11 ]; then
        echo "$filtered_bihet_vcfgz exists and non-empty, newer than $filtered_bi_vcfgz."
        exit
    fi
fi

mkdir -p $tmp_dir

now=$(date +"%T")
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with ${sample} :"${now}

job_ids=()
for N in {1..22} X
do
    log="log/${sample}_step0-${N}-extractVCF-%j.out"

    input_vcfgz=$vcf_dir/ALL.chr${N}.*.vcf.gz
    output_vcf=$tmp_dir/${sample}_chr${N}.vcf

    output_header_file=""
    if [ $N == 21 ]; then
        output_header_file=$tmp_dir/${sample}_header.vcf
    else
        output_header_file=""
    fi

    job_id=$(sbatch \
        --parsable \
        --job-name=${sample}-extractVCF-chr${N} \
        --out=$log --error=$log \
        $scripts_dir/0__extract_chrom.slurm \
        $sample $input_vcfgz $output_vcf $output_header_file \
    )
    job_ids+=( $job_id )
done

log="log/${sample}_combineVCF-%j.out"
combine_job_id=$(sbatch \
    --parsable \
    --kill-on-invalid-dep=yes \
    --dependency=afterok:$(echo "${job_ids[@]}" | tr ' ' ':') \
    --job-name=${SLURM_JOB_NAME}-combine \
    --out=$log --error=$log \
    $scripts_dir/0__combine_VCF.slurm $sample $tmp_dir $bi_vcfgz $bihet_vcfgz \
)

$scripts_dir/wait_for_jobs.sh ${job_ids[@]} "$combine_job_id" || exit 1

now=$(date +"%T")
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> done with ${sample} :"${now}
