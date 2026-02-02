import tempfile
import argparse
from cyvcf2 import VCF, Writer
import os
import glob
from tabulate import tabulate


# Default chromosomes (1-22)
DEFAULT_CHROMS = ["chr"+str(i) for i in range(1, 23)] + ['chrX']

####################
# for debugging
def cython_obj_to_dict(obj: object) -> dict:
    keys = dir(obj)
    keys = filter(lambda k: k[0] != "_" and not k.isupper(), keys)
    data = [(k, getattr(obj, k, None)) for k in keys]
    data = [(k, v) for k, v in data if not callable(v)]
    data = dict(data)
    return data

def display_obj_values(obj: object):
    data = cython_obj_to_dict(obj)
    print(tabulate(data.items()))
####################

def strip_chr(chrom):
    """Remove 'chr' prefix from chromosome name if present."""
    return chrom[3:] if str(chrom).lower().startswith('chr') else str(chrom)


# vcf = vcf object from cyvcf2
def filter_by_chrom_and_sort(vcf, chroms=DEFAULT_CHROMS):
    # Create a temporary file
    temp_fd, temp_path = tempfile.mkstemp(suffix='.vcf')
    os.close(temp_fd)
    w = Writer(temp_path, vcf)
    # Filter and sort
    chrom_order = {chrom: idx for idx, chrom in enumerate(chroms)}
    print("filter chrom!")
    records = [variant for variant in vcf if variant.CHROM in set(chroms)]
    print("sort!")
    records.sort(key=lambda v: (chrom_order[v.CHROM], v.POS))
    # Write sorted records to temp file
    for record in records:
        w.write_record(record)
    w.close()
    return VCF(temp_path)


# write_mode "w": uncompressed VCF, "wz" compressed VCF
def filter_vcf(vcf_in, vcf_out, sv_size_min=40, write_mode='w', remove_chr=False):
    vcf = VCF(vcf_in)
    vcf = filter_by_chrom_and_sort(vcf)

    w = Writer(vcf_out, vcf, "w")  

    print("filter size!")
    for variant in vcf:
        # Get variant type
        sv_type = variant.INFO.get('SVTYPE', '')

        # Initialize pass_filters flag
        pass_filters = True
        if sv_type == 'SNV':
            continue

        # filter imprecise (Sniffles2-specific)
        is_sniffles = 'sniffles' in variant.ID.lower()
        if is_sniffles:
            if not variant.INFO.get('PRECISE', False):
                pass_filters = False

        # Filter by size
        if sv_type != "SNV":
            # Get SV size from SVLEN
            sv_size = abs(variant.INFO.get('SVLEN', 0))
            if isinstance(sv_size, tuple):
                sv_size = abs(sv_size[0])
            if sv_size < sv_size_min:
                pass_filters = False

        # Write variant if it passes filters
        if pass_filters:
            if remove_chr:
                # Remove 'chr' from chromosome name before writing
                variant.CHROM = strip_chr(variant.CHROM)
            w.write_record(variant)

    w.close()
    vcf.close()


def main():
    parser = argparse.ArgumentParser(description='Filter VCFs by chromosome and length')
    parser.add_argument('vcf_files', help='Input VCF file(s) (supports wildcards)', nargs='+')
    parser.add_argument('--sv-size-min', type=int, default=40,
                      help='Minimum SV size to include (default: 40)')
    parser.add_argument('--write-mode', type=str, default="w", help='cyvcf2 VCF Writer mode (default: w)')
    args = parser.parse_args()

    print("filter for variants in chroms:")
    print(DEFAULT_CHROMS)

    print("minimum SV size to keep: {}".format(args.sv_size_min))

    # Process all input files that match the wildcard pattern
    for pattern in args.vcf_files:
        for vcf_file in glob.glob(pattern):
            print("processing {}".format(vcf_file))
            out_vcf = vcf_file.replace("vcf", "filtered.vcf")
            if not os.path.isfile(out_vcf):
                filter_vcf(vcf_file, out_vcf, args.sv_size_min, args.write_mode)
            else:
                print("output already exists!")


if __name__ == '__main__':
    main()
