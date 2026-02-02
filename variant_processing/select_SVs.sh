#!/bin/bash

set -e


# Convert ref allele sequence to N for matching SVs to VAPOR result
#bcftools view kanpig/kanpig.vcf.gz | \
#    awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$4="N"; print}' | \
#    bcftools view -Oz -o kanpig/kanpig.refN.vcf.gz -
#tabix -p vcf kanpig/kanpig.refN.vcf.gz
#
# Select INV SVs from VAPOR result
#bcftools filter -i 'INFO/SVTYPE == "INV"' vapor/combined_vapor_GT_INV.vcf.gz -Oz -o vapor/vapor.onlyINV.vcf.gz
#tabix -p vcf vapor/vapor.onlyINV.vcf.gz
#
# Update INV GT calls using VAPOR
#bcftools annotate \
#	-a vapor/vapor.onlyINV.vcf.gz \
#	-c INFO/VaPoR_QS,INFO/VaPoR_GS,INFO/VaPoR_Rec,FORMAT/GT,FORMAT/GQ \
#	-Oz \
#	-o updateGT.sv.vcf.gz \
#	kanpig/kanpig.refN.vcf.gz
#
#tabix -p vcf updateGT.sv.vcf.gz



# Use caller provided GT calls for INV

# filter out INV GT calls from Kanpig (not reliable)
bcftools filter -e 'INFO/SVTYPE == "INV"' kanpig/kanpig.vcf.gz -Oz -o kanpig/kanpig.noINV.vcf.gz
tabix -p vcf kanpig/kanpig.noINV.vcf.gz

# take INV GT calls from caller
bcftools filter -i 'INFO/SVTYPE == "INV"' ../add_QC/final.cohort.vcf.gz -Oz -o ../add_QC/final.cohort.onlyINV.vcf.gz
tabix -p vcf ../add_QC/final.cohort.onlyINV.vcf.gz

bcftools concat -a -Oz -o final.cohort.vcf.gz kanpig/kanpig.noINV.vcf.gz ../add_QC/final.cohort.onlyINV.vcf.gz
tabix -p vcf final.cohort.vcf.gz

# Exclude SVs based on GT calling accuracy
#   Acceptable sizes
#     - INS & DEL: 50-10k
#     - DUP: 50-1k
#     - INV: Any
bcftools filter \
	-e 'INFO/SVTYPE != "INV" && ( (ABS(INFO/SVLEN) > 10000) || (INFO/SVTYPE == "DUP" && ABS(INFO/SVLEN) > 1000) )' \
	-Oz \
	-o cohort.size_filtered.vcf.gz \
	final.cohort.vcf.gz
tabix -p vcf cohort.size_filtered.vcf.gz


# Select SVs around 100kb from PD GWAS eGenes
bedtools intersect -a cohort.size_filtered.vcf.gz \
        -b /stanley/levin_asap_storage/kwanho/Revio_gDNA/refs/no_alt_GRCh38_ref/pd_124_genes_100kb_padding.bed \
        -header -u > pd.sv.vcf
bcftools sort pd.sv.vcf -Oz -o pd.sv.sort.vcf.gz
tabix -p vcf pd.sv.sort.vcf.gz


# update INFO tags
bcftools +fill-tags pd.sv.sort.vcf.gz -Oz -o pd.sv.filtered.vcf.gz -- -S table_sampleID_Condition.tsv
tabix -p vcf pd.sv.filtered.vcf.gz


# exclude monomorphic SVs, missing GT calls > 95%, and variants only found in HC
bcftools filter -i 'F_MISSING <= 0.95 && AF < .99 && (AC_PD > 1 || AC_ILBD > 1)' pd.sv.filtered.vcf.gz -O z -o final.pd.sv.vcf.gz
tabix -p vcf final.pd.sv.vcf.gz

