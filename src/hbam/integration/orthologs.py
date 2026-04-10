"""Cross-species ortholog mapping (human <-> mouse)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from hbam.utils.logging import log_filter, log_metric, log_step


# Built-in common ortholog mappings (fallback when BioMart unavailable)
_COMMON_ORTHOLOGS = {
    # Mouse -> Human (common muscle/aging genes)
    "Ttn": "TTN", "Myh1": "MYH1", "Myh2": "MYH2", "Myh4": "MYH4",
    "Acta1": "ACTA1", "Des": "DES", "Vim": "VIM", "Ckm": "CKM",
    "Mb": "MB", "Pax7": "PAX7", "Myod1": "MYOD1", "Myf5": "MYF5",
    "Tp53": "TP53", "Mtor": "MTOR", "Cdkn1a": "CDKN1A",
    "Il6": "IL6", "Tnf": "TNF", "Tgfb1": "TGFB1", "Foxo3": "FOXO3",
    "Sirt1": "SIRT1", "Igf1": "IGF1", "Mstn": "MSTN",
    "Col1a1": "COL1A1", "Col3a1": "COL3A1", "Fn1": "FN1", "Lama2": "LAMA2",
    "Sox9": "SOX9", "Acan": "ACAN", "Col2a1": "COL2A1", "Runx2": "RUNX2",
    "Sp7": "SP7", "Gapdh": "GAPDH", "Actb": "ACTB",
}


def load_ortholog_table(
    source: str = "biomart",
    cache_path: Path | None = None,
) -> pd.DataFrame:
    """Load human-mouse ortholog mapping table.

    Args:
        source: "biomart" to query Ensembl, "cache" to load from file, "builtin" for minimal fallback.
        cache_path: Path to cache/load the table.

    Returns:
        DataFrame with columns: human_gene, mouse_gene, orthology_type, sequence_identity.
    """
    with log_step("load_ortholog_table", source=source):
        if cache_path and Path(cache_path).exists() and source != "builtin":
            logger.info(f"Loading cached ortholog table from {cache_path}")
            return pd.read_csv(cache_path, sep="\t")

        if source == "biomart":
            try:
                table = _query_biomart()
                if cache_path:
                    cache_path = Path(cache_path)
                    cache_path.parent.mkdir(parents=True, exist_ok=True)
                    table.to_csv(cache_path, sep="\t", index=False)
                    logger.info(f"Cached ortholog table to {cache_path}")
                return table
            except Exception as e:
                logger.warning(f"BioMart query failed: {e}. Using built-in fallback.")
                return _builtin_table()
        else:
            return _builtin_table()


def _query_biomart() -> pd.DataFrame:
    """Query Ensembl BioMart for human-mouse orthologs."""
    try:
        from pybiomart import Server

        server = Server(host="http://www.ensembl.org")
        dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]

        result = dataset.query(
            attributes=[
                "external_gene_name",
                "mmusculus_homolog_associated_gene_name",
                "mmusculus_homolog_orthology_type",
                "mmusculus_homolog_perc_id",
            ],
            filters={},
        )

        result.columns = ["human_gene", "mouse_gene", "orthology_type", "sequence_identity"]
        result = result.dropna(subset=["human_gene", "mouse_gene"])
        result = result[result["human_gene"].str.len() > 0]
        result = result[result["mouse_gene"].str.len() > 0]

        log_metric("biomart_ortholog_pairs", len(result))
        return result

    except Exception as e:
        raise RuntimeError(f"BioMart query failed: {e}")


def _builtin_table() -> pd.DataFrame:
    """Create a minimal ortholog table from built-in mappings."""
    rows = []
    for mouse, human in _COMMON_ORTHOLOGS.items():
        rows.append({
            "human_gene": human,
            "mouse_gene": mouse,
            "orthology_type": "ortholog_one2one",
            "sequence_identity": 90.0,
        })

    table = pd.DataFrame(rows)
    log_metric("builtin_ortholog_pairs", len(table))
    return table


def map_orthologs(
    genes: list[str],
    source_species: str,
    target_species: str,
    table: pd.DataFrame,
    strategy: str = "best",
) -> dict[str, str]:
    """Map genes from one species to another using ortholog table.

    Args:
        genes: List of gene symbols to map.
        source_species: Source species ("human" or "mouse").
        target_species: Target species ("human" or "mouse").
        table: Ortholog mapping table.
        strategy: "best" (highest sequence identity) or "all" (keep all).

    Returns:
        Dict mapping source gene -> target gene.
    """
    if source_species == target_species:
        return {g: g for g in genes}

    if source_species == "mouse":
        src_col, tgt_col = "mouse_gene", "human_gene"
    else:
        src_col, tgt_col = "human_gene", "mouse_gene"

    gene_set = set(genes)
    matches = table[table[src_col].isin(gene_set)]

    mapping = {}
    one_to_many = 0

    for src_gene in gene_set:
        candidates = matches[matches[src_col] == src_gene]

        if len(candidates) == 0:
            continue
        elif len(candidates) == 1:
            mapping[src_gene] = candidates.iloc[0][tgt_col]
        else:
            one_to_many += 1
            if strategy == "best" and "sequence_identity" in candidates.columns:
                best = candidates.loc[candidates["sequence_identity"].idxmax()]
                mapping[src_gene] = best[tgt_col]
            else:
                mapping[src_gene] = candidates.iloc[0][tgt_col]

    unmapped = len(gene_set) - len(mapping)
    log_metric("ortholog_mapped", len(mapping))
    log_metric("ortholog_unmapped", unmapped)
    log_metric("ortholog_one_to_many", one_to_many)

    return mapping


def harmonize_gene_names(
    adata,
    species: str,
    table: pd.DataFrame,
    target_species: str = "human",
) -> "anndata.AnnData":
    """Convert gene names to a common namespace (default: human HGNC).

    Args:
        adata: Input AnnData.
        species: Species of the input data.
        table: Ortholog mapping table.
        target_species: Target namespace species.

    Returns:
        AnnData with harmonized gene names.
    """
    import anndata as ad

    with log_step("harmonize_gene_names", source=species, target=target_species):
        result = adata.copy()
        original_names = result.var_names.tolist()

        if species == target_species:
            # Just uppercase for consistency
            new_names = [g.upper() for g in original_names]
        else:
            mapping = map_orthologs(original_names, species, target_species, table)
            new_names = []
            for name in original_names:
                if name in mapping:
                    new_names.append(mapping[name].upper())
                else:
                    new_names.append(name.upper())

        result.var["original_gene_name"] = original_names
        result.var["source_species"] = species

        # Make unique
        seen = {}
        unique = []
        for n in new_names:
            if n in seen:
                seen[n] += 1
                unique.append(f"{n}_{seen[n]}")
            else:
                seen[n] = 0
                unique.append(n)

        result.var_names = unique

        log_filter("harmonize", before=len(original_names), after=len(unique), reason="gene name harmonization")

        return result
