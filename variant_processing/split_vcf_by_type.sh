#!/bin/bash

set -e

SV_TYPES=("INS" "DEL" "DUP" "INV")

if ! command -v bcftools &> /dev/null || ! command -v bgzip &> /dev/null || ! command -v tabix &> /dev/null; then
    echo "Error: bcftools, bgzip, and/or tabix are not installed or not in your PATH." >&2
    echo "Please install bcftools and htslib to proceed." >&2
    exit 1
fi

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 <manifest_file> [split_output_directory]" >&2
    echo "  - manifest_file: File with paths to uncompressed VCF files." >&2
    echo "  - split_output_directory (optional): Directory for split files. Defaults to the current directory." >&2
    exit 1
fi

manifest_file="$1"
split_dir_arg="$2" # This might be empty, which is fine.

if [ ! -f "${manifest_file}" ]; then
    echo "Error: Manifest file not found at '${manifest_file}'" >&2
    exit 1
fi

# use provided output directory. If not given, use the current directory.
split_output_dir="${split_dir_arg:-.}"

# Create the target directory for the split files.
# The '-p' flag ensures no error is thrown if the directory already exists.
echo "Split files will be placed in: $(realpath "${split_output_dir}")"
mkdir -p "${split_output_dir}"


while IFS= read -r vcf_file || [[ -n "${vcf_file}" ]]; do
    # Skip empty or blank lines in the manifest file.
    [ -z "${vcf_file}" ] && continue

    # Check that the input file exists before trying to process it.
    if [ ! -f "${vcf_file}" ]; then
        echo "Warning: Input file not found: '${vcf_file}'. Skipping." >&2
        continue
    fi
    # Ensure the input is an uncompressed .vcf file
    if [[ "${vcf_file}" != *.vcf ]]; then
        echo "Warning: Input file '${vcf_file}' is not a .vcf file. Skipping." >&2
        continue
    fi

    echo ""
    echo "--- Processing file: ${vcf_file} ---"

    # Compress and Index
    input_dir=$(dirname "${vcf_file}")
    bgzip -f "${vcf_file}"
    compressed_file="${vcf_file}.gz"
    tabix -p vcf -f "${compressed_file}"

    # Split by SVTYPE
    base_name=$(basename "${compressed_file}" .vcf.gz)

    for sv_type in "${SV_TYPES[@]}"; do
        output_file="${split_output_dir}/${base_name}.${sv_type}.vcf"
        bcftools view -i "INFO/SVTYPE=\"${sv_type}\"" -o "${output_file}" -O v "${compressed_file}"
    done
    echo "--- Finished ${vcf_file} ---"

done < "${manifest_file}"

echo ""
echo "DONE"

