export PATH=~/kwanho/git_repos/SURVIVOR/Debug:$PATH
use .bcftools-1.21

# Filter repeat region
COMBINED_VCF=/stanley/levin_asap_storage/kwanho/Revio_gDNA/variant_calls/merge_SVs/sv.combined.vcf.gz
REPEAT_BED=/stanley/levin_asap_storage/kwanho/Revio_gDNA/refs/no_alt_GRCh38_ref/GRCh38_segdup+gaps+centromere.bed
bcftools filter -T ^$REPEAT_BED $COMBINED_VCF -o sv.filtered.vcf.gz
tabix -p vcf sv.filtered.vcf.gz

bcftools filter -R $REPEAT_BED $COMBINED_VCF -o sv.removed.vcf.gz

# Remove chromosome X
bcftools view -t ^chrX -Oz -o sv.noX.vcf.gz sv.filtered.vcf.gz

# Remove genotype columns
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' sv.noX.vcf.gz -o sv.no_header.vcf
bcftools view -h sv.noX.vcf.gz > header_only.txt

# Manual header cleanup -- keeping chr1-22
cat header_only.txt sv.no_header.vcf > sv.finalist.no_gt.vcf
bgzip sv.finalist.no_gt.vcf
tabix -p vcf sv.finalist.no_gt.vcf.gz

bcftools view \
	-r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
	sv.finalist.no_gt.vcf.gz -Oz -o sv.chr1-22.no_gt.vcf.gz

# Transfer caller predicted genotypes
python transfer_gt.py -v sv.chr1-22.no_gt.vcf.gz -m results_map.txt -o cohort.sv.vcf
bgzip cohort.sv.vcf
tabix -p vcf cohort.sv.vcf.gz

# Make unique IDs
bcftools annotate \
	--set-id '%CHROM\_%POS\_%SVTYPE\_%SVLEN\' \
	cohort.sv.vcf.gz \
	-Oz -o final.preQC.sv.vcf.gz


