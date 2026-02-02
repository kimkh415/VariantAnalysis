export PATH=~/kwanho/git_repos/SURVIVOR/Debug:$PATH
use .bcftools-1.21

SURVIVOR merge filelist_sniffles2_INS.txt 1000 1 1 1 0 50 cohort_merged/sv.sniffles2.INS.vcf
SURVIVOR merge filelist_sniffles2_DEL.txt 1000 1 1 1 0 50 cohort_merged/sv.sniffles2.DEL.vcf
SURVIVOR merge filelist_sniffles2_DUP.txt 1000 1 1 1 0 50 cohort_merged/sv.sniffles2.DUP.vcf
SURVIVOR merge filelist_sniffles2_INV.txt 1000 1 1 1 0 50 cohort_merged/sv.sniffles2.INV.vcf

SURVIVOR merge filelist_pbsv_INS.txt 1000 1 1 1 0 50 cohort_merged/sv.pbsv.INS.vcf
SURVIVOR merge filelist_pbsv_DEL.txt 1000 1 1 1 0 50 cohort_merged/sv.pbsv.DEL.vcf
SURVIVOR merge filelist_pbsv_DUP.txt 1000 1 1 1 0 50 cohort_merged/sv.pbsv.DUP.vcf
SURVIVOR merge filelist_pbsv_INV.txt 1000 1 1 1 0 50 cohort_merged/sv.pbsv.INV.vcf

# cue does not call INS by design
SURVIVOR merge filelist_cue2_INS.txt 1000 1 1 1 0 50 cohort_merged/sv.cue2.INS.vcf
SURVIVOR merge filelist_cue2_DEL.txt 1000 1 1 1 0 50 cohort_merged/sv.cue2.DEL.vcf
SURVIVOR merge filelist_cue2_DUP.txt 1000 1 1 1 0 50 cohort_merged/sv.cue2.DUP.vcf
SURVIVOR merge filelist_cue2_INV.txt 1000 1 1 1 0 50 cohort_merged/sv.cue2.INV.vcf

cd cohort_merged

# make SV IDs unique
bcftools annotate --set-id 'sniffles2.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.sniffles2.INS.vcf sv.sniffles2.INS.vcf
bcftools annotate --set-id 'sniffles2.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.sniffles2.DEL.vcf sv.sniffles2.DEL.vcf
bcftools annotate --set-id 'sniffles2.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.sniffles2.DUP.vcf sv.sniffles2.DUP.vcf
bcftools annotate --set-id 'sniffles2.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.sniffles2.INV.vcf sv.sniffles2.INV.vcf

bcftools annotate --set-id 'pbsv.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.pbsv.INS.vcf sv.pbsv.INS.vcf
bcftools annotate --set-id 'pbsv.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.pbsv.DEL.vcf sv.pbsv.DEL.vcf
bcftools annotate --set-id 'pbsv.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.pbsv.DUP.vcf sv.pbsv.DUP.vcf
bcftools annotate --set-id 'pbsv.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.pbsv.INV.vcf sv.pbsv.INV.vcf

bcftools annotate --set-id 'cue2.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.cue2.INS.vcf sv.cue2.INS.vcf
bcftools annotate --set-id 'cue2.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.cue2.DEL.vcf sv.cue2.DEL.vcf
bcftools annotate --set-id 'cue2.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.cue2.DUP.vcf sv.cue2.DUP.vcf
bcftools annotate --set-id 'cue2.%CHROM\_%POS\_%SVTYPE\.%VKX' -O v -o sv.unique.cue2.INV.vcf sv.cue2.INV.vcf

# merge across callers
ls -1 sv.unique.*.INS.vcf > filelist_INS.txt
ls -1 sv.unique.*.DEL.vcf > filelist_DEL.txt
ls -1 sv.unique.*.DUP.vcf > filelist_DUP.txt
ls -1 sv.unique.*.INV.vcf > filelist_INV.txt

SURVIVOR merge filelist_INS.txt 1000 1 1 1 0 50 ../sv.combined.INS.vcf
SURVIVOR merge filelist_DEL.txt 1000 1 1 1 0 50 ../sv.combined.DEL.vcf
SURVIVOR merge filelist_DUP.txt 1000 1 1 1 0 50 ../sv.combined.DUP.vcf
SURVIVOR merge filelist_INV.txt 1000 1 1 1 0 50 ../sv.combined.INV.vcf


cd ../

# sort and merge across SV types
bcftools sort -Oz -o sv.combined.DEL.vcf.gz sv.combined.DEL.vcf
bcftools sort -Oz -o sv.combined.INS.vcf.gz sv.combined.INS.vcf
bcftools sort -Oz -o sv.combined.DUP.vcf.gz sv.combined.DUP.vcf
bcftools sort -Oz -o sv.combined.INV.vcf.gz sv.combined.INV.vcf
tabix -p vcf sv.combined.DEL.vcf.gz
tabix -p vcf sv.combined.INS.vcf.gz
tabix -p vcf sv.combined.DUP.vcf.gz
tabix -p vcf sv.combined.INV.vcf.gz
ls -1 sv.combined.*.vcf.gz > filelist_svtypes.txt
bcftools concat -f filelist_svtypes.txt -a -O z -o sv.combined.vcf.gz
tabix -p vcf sv.combined.vcf.gz


