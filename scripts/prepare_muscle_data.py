"""Prepare PXD047296 muscle proteomics for the HBAM pipeline."""

import re
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger


def prepare_muscle_proteomics():
    """Convert DIA-NN pg_matrix to pipeline-ready format."""
    input_path = Path("data/raw/PXD047296/MuscleDIANNStep3report.pg_matrix.tsv")
    output_dir = Path("data/processed")
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Reading {input_path}")
    df = pd.read_csv(input_path, sep="\t")

    # Identify metadata vs sample columns
    meta_cols = ["Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"]
    sample_cols = [c for c in df.columns if c not in meta_cols]

    # Extract sample IDs from file paths
    sample_map = {}
    for col in sample_cols:
        m = re.search(r"(\d+)\.mzML", col)
        if m:
            sample_id = int(m.group(1))
            sample_map[col] = f"muscle_S{sample_id:03d}"

    # Gene names: take first gene from semicolon-separated list
    gene_names = []
    for g in df["Genes"].fillna(""):
        if ";" in str(g):
            gene_names.append(str(g).split(";")[0].strip())
        elif g:
            gene_names.append(str(g))
        else:
            gene_names.append("UNKNOWN")

    # Make unique
    seen = {}
    unique_genes = []
    for name in gene_names:
        if name in seen:
            seen[name] += 1
            unique_genes.append(f"{name}_{seen[name]}")
        else:
            seen[name] = 0
            unique_genes.append(name)

    # Build expression matrix (genes as rows, samples as columns)
    expr = df[sample_cols].copy()
    expr.index = unique_genes
    expr.columns = [sample_map.get(c, c) for c in sample_cols]

    # Replace 0 with NaN (DIA-NN uses 0 for not detected)
    expr = expr.replace(0, np.nan)

    # Save expression matrix
    expr_path = output_dir / "muscle_proteomics.tsv"
    expr.to_csv(expr_path, sep="\t")
    logger.info(f"Expression matrix saved: {expr_path} ({expr.shape[0]} genes x {expr.shape[1]} samples)")

    # Create sample metadata
    # Paper design: 4 timepoints, 10 mice each, sequential numbering 281-320
    sample_ids = sorted(sample_map.values())
    metadata_rows = []
    for sid in sample_ids:
        sample_num = int(re.search(r"(\d+)", sid).group(1))

        if sample_num <= 290:
            timepoint = 4
            age_group = "young"
        elif sample_num <= 300:
            timepoint = 8
            age_group = "young"
        elif sample_num <= 310:
            timepoint = 12
            age_group = "old"
        else:
            timepoint = 20
            age_group = "old"

        metadata_rows.append({
            "sample_id": sid,
            "timepoint": timepoint,
            "timepoint_months": f"{timepoint}M",
            "condition": age_group,
            "organism": "mouse",
            "tissue": "skeletal_muscle",
            "dataset": "PXD047296",
        })

    meta_df = pd.DataFrame(metadata_rows)
    meta_path = output_dir / "muscle_metadata.tsv"
    meta_df.to_csv(meta_path, sep="\t", index=False)
    logger.info(f"Metadata saved: {meta_path}")

    # Print summary
    print(f"\n{'='*60}")
    print(f"Muscle Proteomics Data Preparation Summary")
    print(f"{'='*60}")
    print(f"Source: PXD047296 (Ten Organs Proteome Atlas)")
    print(f"Organism: Mus musculus")
    print(f"Tissue: Skeletal muscle")
    print(f"Proteins: {expr.shape[0]}")
    print(f"Samples: {expr.shape[1]}")
    print(f"Timepoints: 4M (n=10), 8M (n=10), 12M (n=10), 20M (n=10)")
    print(f"Conditions: young (4M+8M, n=20), old (12M+20M, n=20)")
    print(f"Missing values: {expr.isna().sum().sum()} ({expr.isna().sum().sum()/expr.size*100:.1f}%)")
    print(f"\nOutput files:")
    print(f"  {expr_path}")
    print(f"  {meta_path}")

    return expr_path, meta_path


if __name__ == "__main__":
    prepare_muscle_proteomics()
