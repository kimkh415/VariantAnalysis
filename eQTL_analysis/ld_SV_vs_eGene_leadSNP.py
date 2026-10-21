# Linkage disequilibrium (r^2 and D') between each significant SV and its eGene lead eSNP from genotype dosages; writes a per-pair LD table.
import argparse
import sys
import numpy as np
import pandas as pd
from cyvcf2 import VCF

UNKNOWN = 3


def strip_chr(c):
    c = str(c)
    return c[3:] if c.startswith("chr") else c


def dosage_map(variant, samples):
    d = variant.gt_types
    return {s: d[i] for i, s in enumerate(samples)}


def ld_stats(d1, d2, common):
    a = np.array([d1[s] for s in common], float)
    b = np.array([d2[s] for s in common], float)
    keep = (a != UNKNOWN) & (b != UNKNOWN)
    a, b = a[keep], b[keep]
    n = len(a)
    if n < 3 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan, np.nan
    r2 = float(np.corrcoef(a, b)[0, 1] ** 2)
    pA, pB = a.sum() / (2 * n), b.sum() / (2 * n)
    nAB_det = 2 * np.sum((a == 2) & (b == 2)) + np.sum((a == 2) & (b == 1)) + np.sum((a == 1) & (b == 2))
    n11 = np.sum((a == 1) & (b == 1))
    fAB = pA * pB
    lo, hi = max(0.0, pA + pB - 1), min(pA, pB)
    for _ in range(100):
        fab, fAb, faB = 1 - pA - pB + fAB, pA - fAB, pB - fAB
        denom = fAB * fab + fAb * faB
        g = (fAB * fab) / denom if denom > 0 else 0.5
        fAB = min(max((nAB_det + 2 * n11 * g) / (2 * n), lo), hi)
    D = fAB - pA * pB
    Dmax = min(pA * (1 - pB), (1 - pA) * pB) if D > 0 else min(pA * pB, (1 - pA) * (1 - pB))
    dprime = min(abs(D) / Dmax, 1.0) if Dmax > 1e-12 else np.nan
    return r2, float(dprime)


def load_lead_lut(lut_tsv):
    df = pd.read_csv(lut_tsv, sep="\t", comment="#")
    df.columns = [c.strip() for c in df.columns]
    return df  # gene, lead_eSNP, chr, bp, source, n_candidates, mega, use


def index_sv_vcf(sv_vcf):
    vcf = VCF(sv_vcf, gts012=True)
    samples = vcf.samples
    sv = {v.ID: (v.CHROM, dosage_map(v, samples)) for v in vcf if v.ID}
    vcf.close()
    return sv, samples


def fetch_snp(snp_vcf, snp_samples, chrom, bp, rsid):
    reg = f"{strip_chr(chrom)}:{int(bp)}-{int(bp)}"
    for v in snp_vcf(reg):
        if v.ID == rsid or int(v.POS) == int(bp):
            return dosage_map(v, snp_samples)
    return None


def main():
    p = argparse.ArgumentParser(description="r^2/D' between each SV and its eGene lead eSNP (single-eSNP lookup).")
    p.add_argument("--pairs", required=True)
    p.add_argument("--lead_snp_lut", required=True, help="lead_esnp_lookup_perGene.tsv")
    p.add_argument("--sv_vcf", required=True)
    p.add_argument("--snp_vcf", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t")
    lut = load_lead_lut(args.lead_snp_lut)
    by_gene = {r.gene: r for r in lut.itertuples(index=False)}

    sv_dos, sv_samples = index_sv_vcf(args.sv_vcf)
    snp_vcf = VCF(args.snp_vcf, gts012=True)
    snp_samples = snp_vcf.samples
    common = [s for s in sv_samples if s in set(snp_samples)]
    if len(common) < 3:
        sys.exit(f"Too few shared samples (n={len(common)}).")

    rows = []
    for row in pairs.itertuples(index=False):
        sv_id, gene = row.SV_id, row.Gene
        rec = by_gene.get(gene)
        if rec is None or sv_id not in sv_dos:
            rows.append((sv_id, gene, None, np.nan, np.nan, "no_lead_snp_or_sv"))
            continue
        if str(rec.use).strip().lower() != "true":       # eSNP absent from MEGA -> excluded, no proxy
            rows.append((sv_id, gene, rec.lead_eSNP, np.nan, np.nan, "lead_snp_excluded"))
            continue
        chrom, sv_d = sv_dos[sv_id]
        snp_d = fetch_snp(snp_vcf, snp_samples, rec.chr, rec.bp, rec.lead_eSNP)
        if snp_d is None:
            rows.append((sv_id, gene, rec.lead_eSNP, np.nan, np.nan, "snp_not_found"))
            continue
        r2v, dpv = ld_stats(sv_d, snp_d, common)
        rows.append((sv_id, gene, rec.lead_eSNP, r2v, dpv, "ok"))
    snp_vcf.close()

    out = pd.DataFrame(rows, columns=["SV_id", "Gene", "lead_GWAS_SNP", "r2", "Dprime", "status"])
    out.to_csv(args.out, sep="\t", index=False)
    vc = out["status"].value_counts().to_dict()
    print(f"Wrote {args.out} ({len(out)} pairs); status: {vc}")


if __name__ == "__main__":
    main()
