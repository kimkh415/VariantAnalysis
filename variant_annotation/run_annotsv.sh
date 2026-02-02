#! /bin/bash

#$ -cwd

#$ -q broad
#$ -P regevlab
#$ -l h_vmem=16g
#$ -l h_rt=5:00:00
#$ -l os=RedHat7
#$ -pe smp 2
#$ -binding linear:2
#$ -R y

source /broad/software/scripts/useuse
source /stanley/levin_dr/kwanho/miniforge3/etc/profile.d/conda.sh
export PATH=$HOME/kwanho/miniforge3/condabin:$PATH

eval "$(mamba shell hook --shell bash)"
mamba activate /stanley/levin_dr/kwanho/anaconda3/envs/truvari-env

use UGER
use .bcftools-1.21
use .bedtools-2.29.0

/stanley/levin_dr/kwanho/git_repos/AnnotSV/bin/AnnotSV \
        -SVinputFile ../PD_eGenes/final.cohort.vcf.gz \
        -outputFile ./cohort.AnnotSV.annotated.tsv \
        -REreport 1 \
        -missingGTinSamplesid 0 \
        -genomeBuild GRCh38 \
        -vcf 1

