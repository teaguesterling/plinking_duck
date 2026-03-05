#!/usr/bin/env python3
"""Generate large WES-like test dataset: ~30K variants x 10K samples on chr10.

Simulates whole exome sequencing data with realistic MAF distribution:
  - ~80% rare variants (MAF < 1%)
  - ~15% low-frequency (MAF 1-5%)
  - ~5% common (MAF > 5%)

Requires: numpy, plink2 (in PATH or PLINK2 env var)

Output: test/data/wes_chr10.{pgen,pvar,psam}
"""

import numpy as np
import subprocess
import os
import sys

SEED = 42
N_SAMPLES = 10_000
N_VARIANTS = 30_000
CHROM = "10"
CHR_LEN = 133_797_422  # chr10 GRCh38

BASES = ["A", "C", "G", "T"]
GT_TABLE = np.array([b"\t0/0", b"\t0/1", b"\t1/1"], dtype="|S4")


def generate_mafs(n, rng):
    """Draw MAFs from a WES-like site frequency spectrum."""
    n_rare = int(n * 0.80)
    n_lowfreq = int(n * 0.15)
    n_common = n - n_rare - n_lowfreq

    rare = rng.beta(0.15, 30, n_rare)           # heavily skewed toward 0
    lowfreq = rng.uniform(0.01, 0.05, n_lowfreq)
    common = rng.uniform(0.05, 0.45, n_common)

    mafs = np.concatenate([rare, lowfreq, common])
    rng.shuffle(mafs)

    # Minimum MAF: 1 alt allele in 2*N diploid alleles
    min_maf = 1.0 / (2 * N_SAMPLES)
    mafs = np.maximum(mafs, min_maf)
    return mafs


def generate_positions(n, rng):
    """Generate positions clustered in exon-like regions across chr10."""
    # ~500 exons of 100-500bp spread across the chromosome
    n_exons = 500
    margin = 1_000_000
    exon_starts = np.sort(rng.integers(margin, CHR_LEN - margin, n_exons))
    exon_lengths = rng.integers(100, 500, n_exons)

    # Distribute variants across exons
    per_exon = np.full(n_exons, n // n_exons)
    per_exon[: n - per_exon.sum()] += 1

    positions = []
    for i in range(n_exons):
        offsets = np.sort(rng.integers(0, max(exon_lengths[i], 1), per_exon[i]))
        positions.append(exon_starts[i] + offsets)

    positions = np.unique(np.concatenate(positions))

    # Fill any deficit from deduplication
    while len(positions) < n:
        extra = rng.integers(margin, CHR_LEN - margin, n - len(positions))
        positions = np.unique(np.concatenate([positions, extra]))

    return np.sort(positions[:n])


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    vcf_path = os.path.join(script_dir, "wes_chr10.vcf")
    out_prefix = os.path.join(script_dir, "wes_chr10")
    plink2 = os.environ.get("PLINK2", "plink2")

    rng = np.random.default_rng(SEED)

    # --- Generate MAFs and genotypes ---
    print(f"Generating {N_VARIANTS} variants x {N_SAMPLES} samples...")

    # Over-generate to account for monomorphic sites after sampling
    n_generate = int(N_VARIANTS * 1.2)
    mafs = generate_mafs(n_generate, rng)

    print("  Sampling genotypes from MAF distribution...")
    genos = rng.binomial(2, mafs[:, np.newaxis],
                         (n_generate, N_SAMPLES)).astype(np.uint8)

    # Filter monomorphic sites
    row_sums = genos.sum(axis=1)
    non_mono = (row_sums > 0) & (row_sums < 2 * N_SAMPLES)
    genos = genos[non_mono][:N_VARIANTS]
    mafs = mafs[non_mono][:N_VARIANTS]
    actual_n = len(genos)
    print(f"  {actual_n} non-monomorphic variants retained")

    # --- Positions and alleles ---
    positions = generate_positions(actual_n, rng)
    ref_idx = rng.integers(0, 4, actual_n)
    alt_idx = (ref_idx + rng.integers(1, 4, actual_n)) % 4

    # --- Write VCF ---
    print(f"Writing VCF to {vcf_path}...")
    sample_names = [f"S{i:05d}" for i in range(1, N_SAMPLES + 1)]

    with open(vcf_path, "wb") as f:
        # Header
        header = "##fileformat=VCFv4.3\n"
        header += f"##contig=<ID={CHROM},length={CHR_LEN}>\n"
        header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        for s in sample_names:
            header += f"\t{s}"
        header += "\n"
        f.write(header.encode())

        # Data lines — use numpy for fast genotype encoding
        for i in range(actual_n):
            prefix = (f"{CHROM}\t{positions[i]}\tv{i+1}\t"
                      f"{BASES[ref_idx[i]]}\t{BASES[alt_idx[i]]}\t.\t.\t.\tGT")
            gt_encoded = GT_TABLE[genos[i]]  # (N_SAMPLES,) array of |S4
            f.write(prefix.encode() + gt_encoded.tobytes() + b"\n")

            if (i + 1) % 5000 == 0:
                print(f"  {i+1}/{actual_n} variants written")

    vcf_size_mb = os.path.getsize(vcf_path) / (1024 * 1024)
    print(f"  VCF size: {vcf_size_mb:.0f} MB")

    # --- Convert to pgen ---
    print("Converting to pgen with plink2...")
    result = subprocess.run(
        [plink2, "--vcf", vcf_path, "--make-pgen", "--out", out_prefix],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("plink2 stderr:", result.stderr, file=sys.stderr)
        sys.exit(1)

    os.unlink(vcf_path)

    # --- Report ---
    for ext in [".pgen", ".pvar", ".psam"]:
        path = out_prefix + ext
        size = os.path.getsize(path)
        if size > 1024 * 1024:
            print(f"  {os.path.basename(path)}: {size / (1024*1024):.1f} MB")
        else:
            print(f"  {os.path.basename(path)}: {size / 1024:.0f} KB")

    # MAF summary
    actual_mafs = genos.sum(axis=1) / (2 * N_SAMPLES)
    actual_mafs = np.minimum(actual_mafs, 1 - actual_mafs)
    print(f"\nMAF distribution:")
    print(f"  Singletons (MAC=1):  {(actual_mafs <= 1/(2*N_SAMPLES)).mean():.1%}")
    print(f"  Rare (MAF < 1%):     {(actual_mafs < 0.01).mean():.1%}")
    print(f"  Low-freq (1-5%):     {((actual_mafs >= 0.01) & (actual_mafs < 0.05)).mean():.1%}")
    print(f"  Common (MAF >= 5%):  {(actual_mafs >= 0.05).mean():.1%}")
    print(f"  Median MAF:          {np.median(actual_mafs):.4f}")
    print(f"  Mean MAF:            {np.mean(actual_mafs):.4f}")


if __name__ == "__main__":
    main()
