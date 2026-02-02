for f in *.vcf; do sed -i 's/NA/./g' $f; done
for f in *.vcf; do bgzip $f; done
for f in *.vcf.gz; do tabix -p vcf $f; done

vcf_list="*.vcf.gz"

bcftools merge \
    -m both \
    -o combined_vapor_GT_INV.vcf.gz \
    -Oz \
    ${vcf_list}

tabix -p vcf combined_vapor_GT_INV.vcf.gz

truvari vcf2df -f -i combined_vapor_GT_INV.vcf.gz combined_vapor_GT_INV.jl
