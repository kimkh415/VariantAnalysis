ls -d ${PWD}/BN*.gz > filelist_kanpig.txt
bcftools merge -m none -l filelist_kanpig.txt | bcftools +fill-tags -O z -o kanpig.vcf.gz
tabix -p vcf kanpig.vcf.gz
