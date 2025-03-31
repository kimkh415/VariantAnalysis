version 1.0

workflow LongcallD_SV_Caller {
    input {
        File ref_fa
        File bam
        File longcallD_tarball
        String sample_name = "sample"
        Int threads = 8
        Int memory_gb = 16
    }

    call RunLongcallD {
        input:
            ref_fa = ref_fa,
            bam = bam,
            longcallD_tarball = longcallD_tarball,
            sample_name = sample_name,
            threads = threads,
            memory_gb = memory_gb
    }

    output {
        File sv_vcf = RunLongcallD.sv_vcf
        File sv_vcf_index = RunLongcallD.sv_vcf_index
        File phased_bam = RunLongcallD.phased_bam
        File phased_bai = RunLongcallD.phased_bai
    }
}

task RunLongcallD {
    input {
        File ref_fa
        File bam
        File longcallD_tarball
        String sample_name
        Int threads
        Int memory_gb
    }

    Int disk_size = ceil(size(ref_fa, "GB") + size(bam, "GB") * 2 + 20)

    command <<<
        set -e

        # Install required tools that might be missing
        apt-get update -y
        apt-get install -y build-essential gcc make cmake libcurl4-openssl-dev zlib1g-dev libbz2-dev liblzma-dev libssl-dev samtools tabix

        # Extract the source tarball
        mkdir -p longcallD_source
        tar -xzf ~{longcallD_source_tarball} -C longcallD_source --strip-components=1
        cd longcallD_source

        # Compile longcallD
        mkdir -p build
        cd build
        cmake ..
        make -j~{threads}

        cd ../../

        # Run longcallD with phasing
        echo "Starting LongcallD..."
        ./longcallD_source/build/longcallD call \
            -t ~{threads} \
            ~{ref_fa} \
            ~{bam} \
            --hifi \
            -b ~{sample_name}.phased.bam > ~{sample_name}.longcallD.vcf

        # Index the VCF and phased BAM
        bgzip -f ~{sample_name}.longcallD.vcf
        tabix -p vcf ~{sample_name}.longcallD.vcf.gz
        samtools index ~{sample_name}.phased.bam
    >>>

    output {
        File sv_vcf = "~{sample_name}.longcallD.vcf.gz"
        File sv_vcf_index = "~{sample_name}.longcallD.vcf.gz.tbi"
        File phased_bam = "~{sample_name}.phased.bam"
        File phased_bai = "~{sample_name}.phased.bam.bai"
    }

    runtime {
        docker: "ubuntu:22.04"
        cpu: threads
        memory: memory_gb
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 3
    }
}

