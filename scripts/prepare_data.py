"""Inspect and prepare downloaded data for the HBAM pipeline."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from loguru import logger


DATA_DIR = Path("data/raw")


def inspect_pride_dataset(accession: str) -> dict:
    """Inspect downloaded PRIDE dataset files."""
    dataset_dir = DATA_DIR / accession

    if not dataset_dir.exists():
        logger.warning(f"Dataset directory not found: {dataset_dir}")
        return {}

    files = list(dataset_dir.iterdir())
    logger.info(f"\n{'='*50}")
    logger.info(f"Dataset: {accession}")
    logger.info(f"Files: {len(files)}")

    info = {"accession": accession, "files": [], "usable": []}

    for f in sorted(files):
        size_mb = f.stat().st_size / 1e6
        logger.info(f"  {f.name} ({size_mb:.1f} MB)")
        info["files"].append({"name": f.name, "size_mb": size_mb, "path": str(f)})

        # Try to read tabular files
        if f.suffix in (".txt", ".tsv", ".csv", ".xlsx"):
            try:
                if f.suffix == ".xlsx":
                    df = pd.read_excel(f, nrows=5)
                else:
                    sep = "\t" if f.suffix in (".txt", ".tsv") else ","
                    df = pd.read_csv(f, sep=sep, nrows=5, low_memory=False)

                logger.info(f"    Columns ({len(df.columns)}): {list(df.columns[:10])}")
                logger.info(f"    Shape preview: {df.shape}")

                # Check for proteomics markers
                cols_lower = [c.lower() for c in df.columns]
                is_maxquant = any("intensity" in c for c in cols_lower)
                is_diann = any("pg.maxlfq" in c or "protein.group" in c for c in cols_lower)
                is_matrix = all(df[c].dtype in ["float64", "int64"] for c in df.columns[1:] if c in df.columns)

                if is_maxquant:
                    info["usable"].append({"path": str(f), "format": "maxquant", "type": "proteomics"})
                    logger.info(f"    -> MaxQuant format detected!")
                elif is_diann:
                    info["usable"].append({"path": str(f), "format": "diann", "type": "proteomics"})
                    logger.info(f"    -> DIA-NN format detected!")
                elif is_matrix:
                    info["usable"].append({"path": str(f), "format": "csv", "type": "generic_matrix"})
                    logger.info(f"    -> Generic matrix format")

            except Exception as e:
                logger.debug(f"    Could not read: {e}")

    return info


def inspect_gem_files() -> list[dict]:
    """Find and inspect Stereo-seq GEM files."""
    gem_files = list(DATA_DIR.rglob("*.gem")) + list(DATA_DIR.rglob("*.gem.gz"))

    results = []
    for f in gem_files:
        logger.info(f"\nGEM file: {f}")
        try:
            import gzip
            opener = gzip.open if f.suffix == ".gz" else open
            with opener(f, "rt") as fh:
                header_lines = []
                data_lines = []
                for i, line in enumerate(fh):
                    if line.startswith("#"):
                        header_lines.append(line.strip())
                    else:
                        data_lines.append(line.strip())
                    if len(data_lines) >= 5:
                        break

                if data_lines:
                    cols = data_lines[0].split("\t")
                    logger.info(f"  Columns: {cols}")
                    logger.info(f"  Data rows preview: {len(data_lines) - 1}")
                    results.append({"path": str(f), "format": "gem", "columns": cols})
        except Exception as e:
            logger.warning(f"  Could not inspect: {e}")

    return results


def generate_config_snippet(datasets: list[dict]) -> str:
    """Generate YAML config snippet for discovered datasets."""
    lines = ["data:", "  output_dir: results", "  intermediate_dir: results/intermediate", "  modalities:"]

    for ds in datasets:
        path = ds.get("path", "")
        fmt = ds.get("format", "csv")
        name = Path(path).stem.replace(".", "_").replace(" ", "_")[:30]

        lines.append(f"    - name: {name}")
        lines.append(f'      path: "{path}"')
        lines.append(f"      format: {fmt}")
        lines.append(f"      species: mouse  # VERIFY")
        lines.append(f"      modality_type: temporal  # VERIFY: temporal, spatial, or functional")
        lines.append(f"      condition_col: condition  # VERIFY column name")
        lines.append("")

    return "\n".join(lines)


def main():
    logger.remove()
    logger.add(sys.stderr, level="INFO", format="<green>{time:HH:mm:ss}</green> | {message}")

    print("=" * 60)
    print("HBAM Pipeline - Data Inspection")
    print("=" * 60)

    all_datasets = []

    # Inspect PRIDE datasets
    for accession in ["PXD047296", "PXD011967", "PXD040722"]:
        info = inspect_pride_dataset(accession)
        if info.get("usable"):
            all_datasets.extend(info["usable"])

    # Inspect GEM files
    gem_results = inspect_gem_files()
    all_datasets.extend(gem_results)

    # Generate config
    if all_datasets:
        print("\n" + "=" * 60)
        print("Suggested config snippet:")
        print("=" * 60)
        print(generate_config_snippet(all_datasets))
    else:
        print("\nNo usable datasets found. Run download_data.py first.")
        print("  python scripts/download_data.py")


if __name__ == "__main__":
    main()
