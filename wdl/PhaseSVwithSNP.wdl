version 1.0


workflow PhaseLongRead {
    input{
        File snpVCF
        File svVCF
        File bam
        String vcf_col
        String vcf_col_sv
        File outdir
        String outdirname = "results"
        File jarFilURL="gs://fc-746adc43-182a-48d6-9cf2-602d748d47d8/TestPhase/AlleleMASSeq.jar"
        File pythonScriptURL="gs://fc-746adc43-182a-48d6-9cf2-602d748d47d8/TestPhase/combineVCF.py"
    }


    call PrepVCF_GZ as PrepVCF_SNP
    {
        input:
            vcf=snpVCF,
            vcf_col=vcf_col
    }

    call PrepVCF_GZ as PrepVCF_SV
    {
        input:
            vcf=svVCF,
            vcf_col=vcf_col_sv

    }

    call PhasedReads
    {
        input:
            bam=bam,
            vcf=PrepVCF_SNP.newVCF,
            jarFil=jarFilURL
    }

    call PrepBAM as PrepBAMOne
    {
        input:
            bam=PhasedReads.all1Bam
    }

    call PrepBAM as PrepBAMTwo
    {
        input:
            bam=PhasedReads.all2Bam
    }


    call CallAlleleSV as CallAlleleOneSV
    {
        input:
            bam=PrepBAMOne.filtBam,
            bai=PrepBAMOne.filtBamBai,
            vcf=PrepVCF_SV.newVCF
    }

    call CallAlleleSV as CallAlleleTwoSV
    {
        input:
            bam=PrepBAMTwo.filtBam,
            bai=PrepBAMTwo.filtBamBai,
            vcf=PrepVCF_SV.newVCF
    }

    call CombineVCF 
    {
        input:
            vcf1=CallAlleleOneSV.allele_vcf,
            vcf2=CallAlleleTwoSV.allele_vcf,
            pythonScript=pythonScriptURL
    }

    output{
        File combinedVCF=CombineVCF.phased_vcf
        File Allele1VCF=CallAlleleOneSV.allele_vcf
        File Allele2VCF=CallAlleleTwoSV.allele_vcf
    }


}



task PhasedReads{
    input{
        File bam
        File vcf
        File jarFil
    }    

    command{
        java -jar ~{jarFil} ~{bam} ~{vcf} allele1.bam allele2.bam failed.bam
    }

    output{
        File all1Bam="allele1.bam"
        File all2Bam="allele2.bam"
    }

    runtime{
        docker: "amazoncorretto:11"
        zones: "us-central1-b"
        memory: "20G"
        disks: "local-disk 100 HDD"
        cpu: 1
    }
}


task CombineVCF{
    input{
        File vcf1
        File vcf2
        File pythonScript
    }
    
    command{
        pip install pandas
        python ~{pythonScript} ~{vcf1} ~{vcf2} comb.vcf
    }

    output{
        File phased_vcf="comb.vcf"
    }

    runtime{
        docker: "python:3"
        zones: "us-central1-b"
        memory: "20G"
        disks: "local-disk 50 HDD"
        cpu: 1
    }
}



task PrepVCF{
    input{
        File vcf
        String vcf_col
    }
    ##bcftools sort ~{vcf} > input.sort.vcf
    ##bgzip -c input.sort.vcf > input.sort.vcf.gz
    ##tabix -p vcf input.sort.vcf.gz
    
    command{
        mkdir tmp
        bcftools sort ~{vcf} -T ./tmp -Oz -o input.sort.vcf.gz
        echo Index
        bcftools index input.sort.vcf.gz
        echo Sort
        bcftools view -O v -s ~{vcf_col} input.sort.vcf.gz |  bcftools view -O v -g het  > new.sv.vcf
  
    }

    output{
        File newVCF="new.sv.vcf"
    }

    runtime{
        docker: "staphb/bcftools:1.11"
        zones: "us-central1-b"
        memory: "40G"
        disks: "local-disk 50 HDD"
        cpu: 1
    }
}

task PrepVCF_GZ{
    input{
        File vcf
        String vcf_col
    }
    ##bcftools sort ~{vcf} > input.sort.vcf
    ##bgzip -c input.sort.vcf > input.sort.vcf.gz
    ##tabix -p vcf input.sort.vcf.gz
    
    command{
        bcftools index ~{vcf}
        bcftools view -O v -s ~{vcf_col} ~{vcf} |  bcftools view -O v -g het > new.sv.vcf
  
    }

    output{
        File newVCF="new.sv.vcf"
    }

    runtime{
        docker: "staphb/bcftools:1.11"
        zones: "us-central1-b"
        memory: "40G"
        disks: "local-disk 100 HDD"
        cpu: 1
    }
}


task PrepBAM{
    input{
        File bam
    }

    command{
        samtools view -b --expr '[AL]>3 || [AM]>3' ~{bam} > allele1.top.bam
        samtools index allele1.top.bam
    }

    output{
        File filtBam="allele1.top.bam"
        File filtBamBai="allele1.top.bam.bai"
    }

    runtime{
        docker: "staphb/samtools:1.19"
        memory: "50G"
        disks: "local-disk 100 HDD"
        cpu: 1
    }
}


task CallAlleleSV{
    input{
        File bam
        File bai
        File vcf
    }
    
    command{
        sniffles --input ~{bam} --genotype-vcf ~{vcf} --vcf allele_sv.vcf
    }

    output{
        File allele_vcf="allele_sv.vcf"
    }

    runtime{
        docker: "us.gcr.io/broad-dsp-lrma/lr-sniffles2:2.0.6"
        memory: "50G"
        disks: "local-disk 100 HDD"
    }
}
