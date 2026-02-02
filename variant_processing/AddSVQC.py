import argparse
import sys
import numpy as np
import pysam
from cyvcf2 import VCF
from concurrent.futures import ProcessPoolExecutor, as_completed


def remove_chr_prefix(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def calculate_ld_r2(genotypes1, genotypes2):
    """Compute pairwise LD R^2 between two genotype vectors."""
    valid_mask = (genotypes1 != 2) & (genotypes2 != 2)  # exclude UNKNOWN (./.)
    if np.sum(valid_mask) < 3:
        return np.nan
    g1, g2 = genotypes1[valid_mask], genotypes2[valid_mask]
    corr = np.corrcoef(g1, g2)[0, 1]
    return corr**2


def calculate_metrics(
    sv_id, chrom, start, end, gt_types, sample_ids, phenotype_map,
    sample_to_bam_path, snp_vcf_path, flanking_bp, ld_r2_threshold, case_groups, control_groups, min_maf, debug
):
    """Compute variant-level metrics for one structural variant."""
    try:
        # === Genotype interpretation ===
        # cyvcf2.gt_types codes: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
        called_mask = gt_types != 2
        alt_mask = np.isin(gt_types, (1, 3))  # carriers: any ALT copy
        # --- Allele frequency and missingness ---
        maf = np.mean(alt_mask[called_mask]) if np.any(called_mask) else 0.0
        missing_rate = 1.0 - np.mean(called_mask)
        # --- Group-level frequencies ---
        case_indices = [i for i, sid in enumerate(sample_ids)
                        if phenotype_map.get(sid) in case_groups]
        control_indices = [i for i, sid in enumerate(sample_ids)
                           if phenotype_map.get(sid) in control_groups]
        case_rate = np.mean(alt_mask[case_indices]) if case_indices else 0.0
        control_rate = np.mean(alt_mask[control_indices]) if control_indices else 0.0
        differential_rate = case_rate - control_rate
        # --- Identify carriers per group for output ---
        group_variant_carriers = {}
        for i, sid in enumerate(sample_ids):
            group = phenotype_map.get(sid)
            if gt_types[i] in (1, 3) and group:
                group_variant_carriers.setdefault(group, []).append(sid)
        # --- LD metrics ---
        snp_vcf = None
        ld_snps_count, common_snp_count, ld_snp_rate = 0, 0, 0.0
        fetch_start = max(0, start - flanking_bp)
        fetch_end = max(fetch_start, end + flanking_bp)
        try:
            snp_vcf = VCF(snp_vcf_path)
            for snp in snp_vcf(f"{remove_chr_prefix(chrom)}:{fetch_start}-{fetch_end}"):
                snp_gts = snp.gt_types
                if snp_gts is None or np.all(snp_gts == 2):
                    continue
                maf_snp = np.mean(np.isin(snp_gts, (1, 3)))
                if maf_snp < min_maf:
                    continue
                common_snp_count += 1
                r2 = calculate_ld_r2(gt_types, snp_gts)
                if not np.isnan(r2) and r2 >= ld_r2_threshold:
                    ld_snps_count += 1
        except Exception as e:
            if debug:
                print(f"Warning: SNP fetch failed for {chrom}:{fetch_start}-{fetch_end} ({e})", file=sys.stderr)
        finally:
            if snp_vcf is not None:
                snp_vcf.close()
        if common_snp_count > 0:
            ld_snp_rate = ld_snps_count / common_snp_count
        # --- BAM-based read metrics ---
        all_mapping_qualities = []
        reads_per_sample = {}
        variant_carrying_samples = [sample_ids[i] for i in range(len(gt_types)) if gt_types[i] in (1, 3)]
        for sid in variant_carrying_samples:
            bam_path = sample_to_bam_path.get(sid)
            if not bam_path:
                continue
            read_count = 0
            try:
                with pysam.AlignmentFile(bam_path, "rb") as bam_file:
                    for read in bam_file.fetch(chrom, start, end):  # strict interval
                        all_mapping_qualities.append(read.mapping_quality)
                        read_count += 1
            except Exception as e:
                if debug:
                    print(f"Warning: BAM fetch failed for {sid} ({chrom}:{start}-{end}) - {e}", file=sys.stderr)
                continue
        avg_mq = np.mean(all_mapping_qualities) if all_mapping_qualities else 0.0
        total_reads = len(all_mapping_qualities)
        mean_reads_per_sample = np.mean(list(reads_per_sample.values())) if reads_per_sample else 0.0

        return (
            maf,
            missing_rate,
            differential_rate,
            ld_snps_count,
            avg_mq,
            group_variant_carriers,
            case_rate,
            control_rate,
            total_reads,
            mean_reads_per_sample,
            common_snp_count,
            ld_snp_rate,
        )
    except Exception as e:
        print(f"Error processing {sv_id}: {e}", file=sys.stderr)
        return (0, 0, 0, 0, 0, {}, 0, 0, 0, 0, 0)


def main():
    parser = argparse.ArgumentParser(description="Compute structural variant metrics for filtering/scoring.")
    parser.add_argument("--sv_vcf", required=True, help="Input SV VCF file")
    parser.add_argument("--snp_vcf", required=True, help="Input common SNP VCF file")
    parser.add_argument("--bam_map", required=True, help="TSV: sample_id\tpath_to_bam")
    parser.add_argument("--phenotype_map", required=True, help="TSV: sample_id\tphenotype_label")
    parser.add_argument("--case_groups", required=True, help="Comma-separated list of phenotype groups as 'case'")
    parser.add_argument("--control_groups", required=True, help="Comma-separated phenotype groups as 'control'")
    parser.add_argument("--output_vcf", required=True, help="Output annotated VCF")
    parser.add_argument("--min_maf", type=float, default=0.0)
    parser.add_argument("--ld_r2_threshold", type=float, default=0.8)
    parser.add_argument("--flanking_bp", type=int, default=10000)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    # Parse phenotype group arguments
    case_groups = [g.strip() for g in args.case_groups.split(",") if g.strip()]
    control_groups = [g.strip() for g in args.control_groups.split(",") if g.strip()]
    all_groups = case_groups + control_groups

    # Load sample-to-bam map
    sample_to_bam_path = {}
    with open(args.bam_map) as f:
        for line in f:
            if not line.strip():
                continue
            sid, path = line.strip().split("\t")
            sample_to_bam_path[sid] = path

    # Load phenotype map
    phenotype_map = {}
    with open(args.phenotype_map) as f:
        for line in f:
            if not line.strip():
                continue
            sid, group = line.strip().split("\t")
            phenotype_map[sid] = group.strip()

    print("Loading variants...")
    sv_vcf = VCF(args.sv_vcf)
    variants = []
    for sv in sv_vcf:
        variants.append((
            sv.ID or f"{sv.CHROM}_{sv.start}",
            sv.CHROM,
            sv.start,
            sv.INFO.get("END", sv.end),
            sv.gt_types,
            str(sv).strip()
        ))
    samples = sv_vcf.samples
    raw_header_lines = sv_vcf.raw_header.strip().splitlines()
    sv_vcf.close()
    print(f"Loaded {len(variants)} variants. Starting processing...")

    # Process in parallel
    future_to_sv = {}
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        for (sv_id, chrom, start, end, gt_types, sv_str) in variants:
            future = executor.submit(
                calculate_metrics,
                sv_id,
                chrom,
                start,
                end,
                gt_types,
                samples,
                phenotype_map,
                sample_to_bam_path,
                args.snp_vcf,
                args.flanking_bp,
                args.ld_r2_threshold,
                case_groups,
                control_groups,
                args.min_maf,
                args.debug
            )
            future_to_sv[future] = (sv_id, chrom, start, sv_str)

        # Collect results
        with open(args.output_vcf, "w") as out_vcf:
            for line in raw_header_lines:  
                if line.startswith("#CHROM"):  
                     continue  
                out_vcf.write(line + "\n")
            out_vcf.write("##source=add_metric_score_fast_v2.py\n")
            out_vcf.write("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency among called samples\">\n")
            out_vcf.write("##INFO=<ID=MISSING_RATE,Number=1,Type=Float,Description=\"Fraction of missing genotypes across samples\">\n")
            out_vcf.write("##INFO=<ID=DIFFERENTIAL_RATE,Number=1,Type=Float,Description=\"Absolute difference in carrier rate between case and control groups\">\n")
            out_vcf.write("##INFO=<ID=CASE_RATE,Number=1,Type=Float,Description=\"Carrier fraction in case group(s)\">\n")
            out_vcf.write("##INFO=<ID=CONTROL_RATE,Number=1,Type=Float,Description=\"Carrier fraction in control group(s)\">\n")
            out_vcf.write("##INFO=<ID=LD_SNPS_COUNT,Number=1,Type=Integer,Description=\"Number of SNPs in strong LD (r2>=threshold)\">\n")
            out_vcf.write("##INFO=<ID=TOTAL_SNPS_NEARBY,Number=1,Type=Integer,Description=\"Total number of SNPs checked in flanking region\">\n")
            out_vcf.write("##INFO=<ID=LD_SNP_RATE,Number=1,Type=Float,Description=\"Fraction of SNPs with strong LD to SV\">\n")
            out_vcf.write("##INFO=<ID=AVG_MAP_QUALITY,Number=1,Type=Float,Description=\"Average read mapping quality for carriers\">\n")
        out_vcf.write('##INFO=<ID=AVG_READS_VALID_SAMPLE,Number=1,Type=Integer,Description="Average number of reads across samples with valid genotype call">\n')
            out_vcf.write("##INFO=<ID=TOTAL_READS,Number=1,Type=Integer,Description=\"Total number of reads overlapping variant among carriers\">\n")
            for group in all_groups:
                group_upper = group.upper()
                desc = f"Sample IDs carrying the variant in group {group_upper}"
                out_vcf.write(f"##INFO=<ID={group_upper}_CARRIERS,Number=.,Type=String,Description=\"{desc}\">\n")

            # Add column header and samples
            header_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples
            out_vcf.write("\t".join(header_cols) + "\n")

            for future in as_completed(future_to_sv):
                out_id, out_chrom, out_pos, raw_line = future_to_sv[future]
                (
                    maf,
                    mr,
                    dr,
                    ldc,
                    amq,
                    group_carriers,
                    case_rate,
                    control_rate,
                    total_reads,
                    mean_reads_per_sample,
                    total_snps,
                    ld_snp_rate,
                ) = future.result()

                # Build group carrier strings (dynamic)
                carrier_parts = [f"{g.upper()}_CARRIERS={','.join(samp)}" for g, samp in group_carriers.items() if samp]
                new_info_str = (
                    f"MAF={maf:.4f};MISSING_RATE={mr:.4f};DIFFERENTIAL_RATE={dr:.4f};"
                    f"CASE_RATE={case_rate:.4f};CONTROL_RATE={control_rate:.4f};"
                    f"LD_SNPS_COUNT={ldc};TOTAL_SNPS_NEARBY={total_snps};LD_SNP_RATE={ld_snp_rate:.4f};"
                    f"AVG_MAP_QUALITY={amq:.2f};TOTAL_READS={total_reads};"
                    f"AVG_READS_VALID_SAMPLE={mean_reads_per_sample:.2f}"
                    + ("" if not carrier_parts else ";" + ";".join(carrier_parts))
                )
                parts = raw_line.split("\t", 8)
                if len(parts) >= 8:
                    existing_info = parts[7]
                    new_line = "\t".join(parts[:7]) + f"\t{existing_info};{new_info_str}"
                    if len(parts) == 9:
                        new_line += f"\t{parts[8]}"  # append exactly original FORMAT + sample block
                    new_line += "\n"
                    out_vcf.write(new_line)

    print("DONE")


if __name__ == "__main__":
    main()

