version 1.0

workflow KanpigGenotype {
    input {
        # Required inputs
        File input_vcf
        File input_bam
        File reference_fasta
        String sample_name = "SAMPLE"
        String output_prefix = "output"

        # Optional inputs
        File? bed  # A sorted bed file that restricts kanpig to only analyzing variants with starts and ends within a single bed entry
        File? ploidy_bed  # special regions within chromosomes that should have non-diploid genotype
        Int threads = 8

        # Kanpig executable location (GCP path)
        String kanpig_gcp_path
    }

    call KanpigSetup {
        input:
            kanpig_gcp_path = kanpig_gcp_path
    }

    call RunKanpig {
        input:
            kanpig_path = KanpigSetup.kanpig_executable,
            input_vcf = input_vcf,
            input_bam = input_bam,
            reference_fasta = reference_fasta,
            sample_name = sample_name,
            output_prefix = output_prefix,
            ploidy_bed = ploidy_bed,
            threads = threads
    }

    output {
        File output_vcf = RunKanpig.output_vcf
        File output_vcf_index = RunKanpig.output_vcf_index
    }
    output {
         File? KanpigSetup.kanpig_executable
    }
}

task KanpigSetup {
    input {
        String kanpig_gcp_path
        String kanpig_version = "1.1.0"
    }

    command <<<
        set -e

        # Create directory for Kanpig
        mkdir -p /kanpig

        # Check if Kanpig is available at the specified GCP path
        if [[ -n "~{kanpig_gcp_path}" ]]; then
            echo "Checking for Kanpig at ~{kanpig_gcp_path}"

            # Try to copy from GCP path
            if gsutil -q stat "~{kanpig_gcp_path}"; then
                echo "Found Kanpig at ~{kanpig_gcp_path}, copying..."
                gsutil cp "~{kanpig_gcp_path}" /kanpig/kanpig
                chmod +x /kanpig/kanpig
                echo "Using existing Kanpig build"
                exit 0
            else
                echo "No Kanpig found at ~{kanpig_gcp_path}, will install new version"
            fi
        fi

        # If we reach here, we need to install Kanpig
        echo "Setting up Kanpig..."

        # Check if we're installing latest or specific version
        if [[ "~{kanpig_version}" == "latest" ]]; then
            # Get the latest release URL
            RELEASE_URL=$(curl -s https://api.github.com/repos/ACEnglish/kanpig/releases/latest | grep browser_download_url | grep x86_64-unknown-linux | cut -d '"' -f 4)
        else
            RELEASE_URL="https://github.com/ACEnglish/kanpig/releases/download/v~{kanpig_version}/kanpig-v~{kanpig_version}-x86_64-unknown-linux-musl.tar.gz"
        fi

        echo "Downloading from: $RELEASE_URL"

        # Try downloading prebuilt binary first
        if curl -L -s -o /tmp/kanpig.tar.gz "$RELEASE_URL"; then
            echo "Downloaded prebuilt binary"
            tar -xzf /tmp/kanpig.tar.gz -C /tmp
            mv /tmp/kanpig /kanpig/
            chmod +x /kanpig/kanpig
        else
            echo "Failed to download prebuilt binary, building from source"

            # Ensure build dependencies are installed
            apt-get update && apt-get install -y curl build-essential git

            # Install Rust if not already installed
            if ! command -v cargo &> /dev/null; then
                echo "Installing Rust..."
                curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
                source $HOME/.cargo/env
            fi

            # Clone and build Kanpig
            git clone https://github.com/ACEnglish/kanpig.git /tmp/kanpig
            cd /tmp/kanpig
            if [[ "~{kanpig_version}" != "latest" ]]; then
                git checkout v~{kanpig_version} || echo "Version not found, using latest"
            fi
            cargo build --release
            cp target/release/kanpig /kanpig/
        fi

        # Verify installation
        /kanpig/kanpig --version

        # Copy to GCP path if specified and successful
        if [[ -n "~{kanpig_gcp_path}" && -f "/kanpig/kanpig" ]]; then
            echo "Copying Kanpig to ~{kanpig_gcp_path} for future use"
            gsutil cp /kanpig/kanpig "~{kanpig_gcp_path}"
        fi
    >>>

    output {
        String kanpig_executable = "/kanpig/kanpig"
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
        memory: "4 GB"
        disks: "local-disk 10 SSD"
        cpu: 2
    }
}

task RunKanpig {
    input {
        String kanpig_path
        File input_vcf
        File input_bam
        File reference_fasta
        String sample_name
        String output_prefix

        File? input_vcf_index
        File? input_bam_index
        File? reference_fasta_index
        File? ploidy_bed
        Int threads

        # Runtime parameters
        Int memory_gb = 16
        Int disk_size_gb = 100
    }

    String output_vcf_filename = output_prefix + ".vcf.gz"

    command <<<
        set -e

        # Use provided Kanpig path
        export PATH="$(dirname ~{kanpig_path}):$PATH"

        # Create indexes if they don't exist
        if [ ! -f "~{input_vcf}.tbi" ] && [ -z "~{input_vcf_index}" ]; then
            tabix -p vcf ~{input_vcf}
        fi

        if [ ! -f "~{input_bam}.bai" ] && [ -z "~{input_bam_index}" ]; then
            samtools index ~{input_bam}
        fi

        if [ ! -f "~{reference_fasta}.fai" ] && [ -z "~{reference_fasta_index}" ]; then
            samtools faidx ~{reference_fasta}
        fi

        # Run Kanpig genotype
        kanpig gt \
            --input ~{input_vcf} \
            --reads ~{input_bam} \
            --reference ~{reference_fasta} \
            --out ~{output_vcf_filename} \
            --sample ~{sample_name} \
            --threads ~{threads} \
            ~{if defined(ploidy_bed) then "--ploidy-bed " + ploidy_bed else ""}

        # Index the output VCF
        tabix -p vcf ~{output_vcf_filename}
    >>>

    output {
        File output_vcf = output_vcf_filename
        File output_vcf_index = output_vcf_filename + ".tbi"
    }

    runtime {
        docker: "quay.io/biocontainers/tabix:1.11--hdfd78af_0"  # Contains tabix and samtools
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        cpu: threads
    }
}
