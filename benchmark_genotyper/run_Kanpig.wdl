version 1.0

workflow KanpigWorkflow {
    input {
        File variants
        File variants_tbi
        File bam
        File bai
        File reference
        File fai
        String sample
        Int threads = 8
        Int ram_size_gb = 64
    }

    call RunKanpig {
        input:
            variants = variants,
            variants_tbi = variants_tbi,
            bam = bam,
            bai = bai,
            reference = reference,
            fai = fai,
            sample = sample,
            threads = threads,
            ram_size_gb = ram_size_gb
    }

    output {
        File sorted_output = RunKanpig.sorted_output
        File sorted_output_idx = RunKanpig.sorted_output_idx
    }
}

task RunKanpig {
    input {
        File variants
        File variants_tbi
        File reference
        File bam
        File bai
        File fai
        String sample
        Int threads
        Int ram_size_gb
    }
    String outputs = "~{sample}.kanpig.vcf.gz"
    Int disk_size_gb = 200 + ceil(size(reference,"GB")) + 100*ceil(size(variants,"GB")) + 2*ceil(size(bam,"GB"))
    String workdir = "/mnt/disks/cromwell_root/kanpig_genotyped"  # Move to Google Batch changed the cromwelll root path
    command <<<
        set -euxo pipefail
        mkdir -p ~{workdir}
        cd ~{workdir}
        EFFECTIVE_MEM_GB=~{ram_size_gb}
        EFFECTIVE_MEM_GB=$(( ${EFFECTIVE_MEM_GB} - 4 ))
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} -1 ))
        df -h
        echo "thread: ${N_THREADS}"
        echo "diskusage: ~{disk_size_gb}"
        pwd
        export RUST_BACKTRACE="full"
        /usr/bin/time -v /software/kanpig --input ~{variants} \
            --bam ~{bam} \
            --reference ~{reference} \
            --sample ~{sample} \
            --hapsim 0.97 \
            --chunksize 500 \
            --maxpaths 1000 \
            --gpenalty 0.04 \
            --threads ~{threads}  \
            --out ~{workdir}/tmp.vcf && echo "kanpig ok!" || "kanpig failed!"
        bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G -O z ~{workdir}/tmp.vcf > ~{workdir}/~{outputs} && echo "bcfsort ok!" || "bcfsort failed!"
        tabix -p vcf  ~{workdir}/~{outputs} && echo "tabix ok!" || echo "tabix failed!"
        ls -l ~{workdir}
    >>>
    output {
        File sorted_output = "~{workdir}/~{outputs}"
        File sorted_output_idx = "~{workdir}/~{outputs}.tbi"
    }
    runtime {
        docker: "quay.io/zhengxc93/terra-kanpig:latest"  
        cpu: threads
        memory: ram_size_gb + "GB"  
        disks: "local-disk " + disk_size_gb + " SSD"  
    }
}

