#! /bin/bash

#$ -cwd

#$ -q broad
#$ -P regevlab
#$ -l h_vmem=24g
#$ -l h_rt=1:00:00
#$ -l os=RedHat7
#$ -pe smp 1
#$ -binding linear:1
#$ -R y

#$ -t 1-3

source /broad/software/scripts/useuse
source /stanley/levin_dr/kwanho/miniforge3/etc/profile.d/conda.sh
export PATH=$HOME/kwanho/miniforge3/condabin:$PATH

eval "$(mamba shell hook --shell bash)"
mamba activate /stanley/levin_dr/kwanho/anaconda3/envs/truvari-env

SEEDFILE=data.list.txt ##File with parameters
invcf=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}')
outvcf=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $2}')
colvcf=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $3}')

truvari collapse \
	-i ${invcf} \
	-o ${outvcf} \
	-c ${colvcf} \
	--dup-to-ins \
	--reference /stanley/levin_asap_storage/kwanho/Revio_gDNA/refs/no_alt_GRCh38_ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
	--sizemin 40


