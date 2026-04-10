"""Multi-format data loaders for the HBAM pipeline."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from loguru import logger

from hbam.utils.logging import log_filter, log_metric, log_step


def load_maxquant(
    path: Path,
    intensity_prefix: str = "LFQ intensity",
) -> ad.AnnData:
    """Load MaxQuant proteinGroups.txt into AnnData.

    Filters reverse hits, potential contaminants, and only-identified-by-site.

    Args:
        path: Path to proteinGroups.txt.
        intensity_prefix: Prefix for intensity columns.

    Returns:
        AnnData (samples x proteins) with .var["gene_names"] and .layers["raw"].
    """
    with log_step("load_maxquant", path=str(path)):
        df = pd.read_csv(path, sep="\t", low_memory=False)
        n_initial = len(df)

        # Filter contaminants and reverse hits
        for col, label in [
            ("Reverse", "reverse_hits"),
            ("Potential contaminant", "contaminants"),
            ("Only identified by site", "only_by_site"),
        ]:
            if col in df.columns:
                mask = df[col].fillna("").str.strip() != "+"
                before = len(df)
                df = df[mask]
                log_filter(label, before=before, after=len(df), reason=f"removed {col}")

        log_filter("total_filtering", before=n_initial, after=len(df), reason="all contaminant filters")

        # Extract intensity columns
        intensity_cols = [c for c in df.columns if c.startswith(intensity_prefix)]
        if not intensity_cols:
            raise ValueError(
                f"No columns matching prefix '{intensity_prefix}' found. "
                f"Available: {df.columns[:10].tolist()}"
            )

        sample_names = [c.replace(intensity_prefix, "").strip() for c in intensity_cols]

        # Build expression matrix (samples x proteins)
        X = df[intensity_cols].values.T.astype(np.float64)
        X[X == 0] = np.nan  # MaxQuant uses 0 for missing

        # Gene names
        gene_col = "Gene names" if "Gene names" in df.columns else "Gene.names"
        if gene_col not in df.columns:
            gene_col = next(
                (c for c in df.columns if "gene" in c.lower() and "name" in c.lower()), None
            )

        if gene_col and gene_col in df.columns:
            gene_names = df[gene_col].fillna("").astype(str).tolist()
            # Take first gene if multiple (semicolon-separated)
            gene_names = [
                g.split(";")[0].strip() if g else f"UNKNOWN_{i}"
                for i, g in enumerate(gene_names)
            ]
        else:
            gene_names = [f"PROTEIN_{i}" for i in range(len(df))]

        # Make unique
        gene_names = _make_unique(gene_names)

        protein_ids = df["Protein IDs"].tolist() if "Protein IDs" in df.columns else gene_names

        var = pd.DataFrame(
            {
                "gene_names": gene_names,
                "protein_ids": protein_ids,
            },
            index=gene_names,
        )

        obs = pd.DataFrame({"sample_id": sample_names}, index=sample_names)

        adata = ad.AnnData(X=X, obs=obs, var=var, layers={"raw": X.copy()})

        log_metric("loaded_samples", adata.n_obs)
        log_metric("loaded_proteins", adata.n_vars)

        return adata


def load_diann(
    path: Path,
    quantity_col: str = "PG.MaxLFQ",
    protein_col: str = "Protein.Group",
    run_col: str = "Run",
    gene_col: str = "Genes",
) -> ad.AnnData:
    """Load DIA-NN report.tsv (long format) into AnnData.

    Args:
        path: Path to DIA-NN report file.
        quantity_col: Column with quantitative values.
        protein_col: Column with protein group IDs.
        run_col: Column with run/sample names.
        gene_col: Column with gene names.

    Returns:
        AnnData (samples x proteins).
    """
    with log_step("load_diann", path=str(path)):
        df = pd.read_csv(path, sep="\t", low_memory=False)

        for col in [run_col, protein_col, quantity_col]:
            if col not in df.columns:
                raise ValueError(
                    f"Required column '{col}' not found. Available: {df.columns.tolist()}"
                )

        # Pivot to wide format: samples x proteins
        pivot = df.pivot_table(
            values=quantity_col,
            index=run_col,
            columns=protein_col,
            aggfunc="max",  # Take max for duplicate entries
        )

        X = pivot.values.astype(np.float64)
        X[X == 0] = np.nan

        sample_names = pivot.index.tolist()
        protein_names = pivot.columns.tolist()

        # Map protein groups to gene names
        gene_map: dict[str, str] = {}
        if gene_col in df.columns:
            for _, row in df[[protein_col, gene_col]].drop_duplicates().iterrows():
                pg = row[protein_col]
                gene = (
                    str(row[gene_col]).split(";")[0].strip()
                    if pd.notna(row[gene_col])
                    else pg
                )
                gene_map[pg] = gene

        gene_names = [gene_map.get(p, p) for p in protein_names]
        gene_names = _make_unique(gene_names)

        var = pd.DataFrame(
            {
                "gene_names": gene_names,
                "protein_ids": protein_names,
            },
            index=gene_names,
        )

        obs = pd.DataFrame({"sample_id": sample_names}, index=sample_names)

        adata = ad.AnnData(X=X, obs=obs, var=var, layers={"raw": X.copy()})

        log_metric("loaded_samples", adata.n_obs)
        log_metric("loaded_proteins", adata.n_vars)

        return adata


def load_diann_matrix(
    path: Path,
    gene_col: str = "Genes",
) -> ad.AnnData:
    """Load DIA-NN pg_matrix.tsv (wide format: proteins x samples).

    This is the pre-pivoted protein group matrix output by DIA-NN,
    with metadata columns followed by sample intensity columns.

    Args:
        path: Path to pg_matrix.tsv file.
        gene_col: Column containing gene names.

    Returns:
        AnnData (samples x proteins) with gene names preserved.
    """
    import re

    with log_step("load_diann_matrix", path=str(path)):
        df = pd.read_csv(path, sep="\t", low_memory=False)

        # Identify metadata vs sample columns
        # DIA-NN matrix has: Protein.Group, Protein.Ids, Protein.Names, Genes,
        # First.Protein.Description, then sample columns (file paths)
        meta_cols = [c for c in df.columns if not (
            c.endswith(".mzML") or c.endswith(".raw") or c.endswith(".d")
            or c.endswith(".wiff") or "/" in c or "\\" in c
        )]
        sample_cols = [c for c in df.columns if c not in meta_cols]

        if not sample_cols:
            # Fallback: non-string columns are samples
            meta_cols = [c for c in df.columns if df[c].dtype == object]
            sample_cols = [c for c in df.columns if c not in meta_cols]

        # Extract gene names
        if gene_col in df.columns:
            gene_names = df[gene_col].fillna("").astype(str).tolist()
            gene_names = [g.split(";")[0].strip() if g else f"UNKNOWN_{i}"
                         for i, g in enumerate(gene_names)]
        else:
            gene_names = [f"PROTEIN_{i}" for i in range(len(df))]

        gene_names = _make_unique(gene_names)

        # Build matrix (samples x proteins)
        X = df[sample_cols].values.T.astype(np.float64)
        X[X == 0] = np.nan

        # Clean sample names from file paths
        sample_names = []
        for c in sample_cols:
            m = re.search(r"(\d+)\.\w+$", c)
            if m:
                sample_names.append(f"S{m.group(1)}")
            else:
                sample_names.append(c.split("/")[-1].split("\\")[-1])
        sample_names = _make_unique(sample_names)

        protein_ids = df["Protein.Group"].tolist() if "Protein.Group" in df.columns else gene_names

        var = pd.DataFrame({
            "gene_names": gene_names,
            "protein_ids": protein_ids,
        }, index=gene_names)

        obs = pd.DataFrame({"sample_id": sample_names}, index=sample_names)

        adata = ad.AnnData(X=X, obs=obs, var=var, layers={"raw": X.copy()})

        log_metric("loaded_samples", adata.n_obs)
        log_metric("loaded_proteins", adata.n_vars)

        return adata


def load_matrix(
    path: Path,
    sep: str = "auto",
    gene_col: str = "auto",
    transpose: bool | None = None,
) -> ad.AnnData:
    """Load a generic CSV/TSV expression matrix into AnnData.

    Auto-detects separator, gene name column, and orientation.

    Args:
        path: Path to CSV/TSV file.
        sep: Separator ("auto" to detect from extension).
        gene_col: Column with gene names ("auto" to detect first non-numeric column).
        transpose: If True, genes are rows. None = auto-detect.

    Returns:
        AnnData (samples x genes).
    """
    with log_step("load_matrix", path=str(path)):
        path = Path(path)

        if sep == "auto":
            sep = "\t" if path.suffix in (".tsv", ".txt", ".gem") else ","

        df = pd.read_csv(path, sep=sep, low_memory=False)

        # Detect gene column
        detected_gene_col: str | None = None
        if gene_col == "auto":
            for col in df.columns:
                if df[col].dtype == object:
                    detected_gene_col = col
                    break
        elif gene_col in df.columns:
            detected_gene_col = gene_col

        if detected_gene_col and detected_gene_col in df.columns:
            gene_names = df[detected_gene_col].astype(str).tolist()
            df = df.drop(columns=[detected_gene_col])
        else:
            gene_names = [f"GENE_{i}" for i in range(len(df))]

        # Remove any other non-numeric columns
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        df = df[numeric_cols]

        X = df.values.astype(np.float64)
        sample_names = df.columns.tolist()

        # Auto-detect orientation: if more rows than columns, likely genes-as-rows
        if transpose is None:
            transpose = X.shape[0] > X.shape[1]

        if transpose:
            X = X.T
            sample_names, gene_names = gene_names, sample_names
            # Repair shapes if the swap produced a mismatch
            if len(gene_names) != X.shape[1]:
                gene_names = [f"GENE_{i}" for i in range(X.shape[1])]
            if len(sample_names) != X.shape[0]:
                sample_names = [f"SAMPLE_{i}" for i in range(X.shape[0])]

        gene_names = _make_unique([str(g) for g in gene_names])
        sample_names = _make_unique([str(s) for s in sample_names])

        var = pd.DataFrame({"gene_names": gene_names}, index=gene_names)
        obs = pd.DataFrame({"sample_id": sample_names}, index=sample_names)

        adata = ad.AnnData(X=X, obs=obs, var=var, layers={"raw": X.copy()})

        log_metric("loaded_samples", adata.n_obs)
        log_metric("loaded_genes", adata.n_vars)

        return adata


def load_stereo_gem(
    path: Path,
    bin_size: int | None = None,
) -> ad.AnnData:
    """Load Stereo-seq GEM format into AnnData.

    GEM format: tab-separated with columns geneID, x, y, MIDCount (or MIDcount, UMICount).

    Args:
        path: Path to GEM file.
        bin_size: If set, aggregate into spatial bins of this size.

    Returns:
        AnnData with spatial coordinates in .obsm["spatial"].
    """
    with log_step("load_stereo_gem", path=str(path)):
        df = pd.read_csv(path, sep="\t", comment="#")

        # Normalize column names to lowercase for matching
        original_cols = df.columns.tolist()
        df.columns = [c.lower() for c in original_cols]

        gene_col_name = next((c for c in df.columns if "gene" in c), None)
        x_col = next((c for c in df.columns if c == "x"), None)
        y_col = next((c for c in df.columns if c == "y"), None)
        count_col = next(
            (c for c in df.columns if "mid" in c or "umi" in c or "count" in c), None
        )

        if not all([gene_col_name, x_col, y_col, count_col]):
            raise ValueError(
                f"Cannot identify required columns. Found: {df.columns.tolist()}"
            )

        # Bin coordinates if requested
        if bin_size:
            df[x_col] = (df[x_col] // bin_size) * bin_size
            df[y_col] = (df[y_col] // bin_size) * bin_size

        # Create bin IDs from coordinates
        df["bin_id"] = df[x_col].astype(str) + "_" + df[y_col].astype(str)

        # Pivot to matrix: bins x genes
        pivot = df.pivot_table(
            values=count_col,
            index="bin_id",
            columns=gene_col_name,
            aggfunc="sum",
            fill_value=0,
        )

        X = pivot.values.astype(np.float32)
        gene_names = pivot.columns.tolist()
        bin_ids = pivot.index.tolist()

        # Extract coordinates
        coords = np.array(
            [[float(b.split("_")[0]), float(b.split("_")[1])] for b in bin_ids],
            dtype=np.float32,
        )

        var = pd.DataFrame({"gene_names": gene_names}, index=gene_names)
        obs = pd.DataFrame(
            {
                "bin_id": bin_ids,
                "total_counts": X.sum(axis=1).astype(int),
                "n_genes_detected": (X > 0).sum(axis=1).astype(int),
            },
            index=bin_ids,
        )

        adata = ad.AnnData(
            X=X,
            obs=obs,
            var=var,
            obsm={"spatial": coords},
            layers={"raw": X.copy()},
        )

        log_metric("loaded_bins", adata.n_obs)
        log_metric("loaded_genes", adata.n_vars)
        log_metric("total_counts", int(X.sum()))

        return adata


def _make_unique(names: list[str]) -> list[str]:
    """Make a list of names unique by appending suffixes."""
    seen: dict[str, int] = {}
    result = []
    for name in names:
        if not name or name == "nan":
            name = "UNNAMED"
        if name in seen:
            seen[name] += 1
            result.append(f"{name}_{seen[name]}")
        else:
            seen[name] = 0
            result.append(name)
    return result
