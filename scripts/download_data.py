"""Download public aging datasets for the HBAM pipeline."""

from __future__ import annotations

import os
import sys
import json
from pathlib import Path
from urllib.parse import urlparse

import requests
from loguru import logger


DATA_DIR = Path("data/raw")
EXTERNAL_DIR = Path("data/external")


def download_file(url: str, dest: Path, chunk_size: int = 8192) -> Path:
    """Download a file with progress reporting."""
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists():
        logger.info(f"Already exists: {dest}")
        return dest

    logger.info(f"Downloading: {url}")
    logger.info(f"  -> {dest}")

    response = requests.get(url, stream=True, timeout=300)
    response.raise_for_status()

    total = int(response.headers.get("content-length", 0))
    downloaded = 0

    with open(dest, "wb") as f:
        for chunk in response.iter_content(chunk_size=chunk_size):
            f.write(chunk)
            downloaded += len(chunk)
            if total > 0:
                pct = downloaded / total * 100
                print(f"\r  Progress: {downloaded / 1e6:.1f} MB / {total / 1e6:.1f} MB ({pct:.0f}%)", end="", flush=True)

    print()
    logger.info(f"Downloaded: {dest} ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


def download_pride_files(accession: str, file_patterns: list[str] | None = None) -> list[Path]:
    """Download files from PRIDE/ProteomeXchange.

    Args:
        accession: PXD accession number.
        file_patterns: Optional list of filename substrings to filter.

    Returns:
        List of downloaded file paths.
    """
    dest_dir = DATA_DIR / accession
    dest_dir.mkdir(parents=True, exist_ok=True)

    # Get file listing from PRIDE API
    api_url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={accession}"
    logger.info(f"Fetching file list from PRIDE: {accession}")

    try:
        response = requests.get(api_url, timeout=60)
        response.raise_for_status()
        files = response.json()
    except Exception as e:
        logger.error(f"Failed to fetch PRIDE file list: {e}")
        # Fallback: try FTP listing
        logger.info("Trying alternative PRIDE API endpoint...")
        alt_url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}/files"
        try:
            response = requests.get(alt_url, timeout=60)
            response.raise_for_status()
            files = response.json()
        except Exception:
            logger.warning(f"Could not fetch file list for {accession}. Manual download may be needed.")
            return []

    downloaded = []

    for file_info in files:
        # Handle different API response formats
        if isinstance(file_info, dict):
            filename = file_info.get("fileName", file_info.get("name", ""))
            download_url = file_info.get("publicFileLocations", [{}])
            if isinstance(download_url, list) and len(download_url) > 0:
                download_url = download_url[0].get("value", "")
            elif isinstance(download_url, str):
                pass
            else:
                # Try FTP
                download_url = file_info.get("downloadLink", file_info.get("ftpLink", ""))
        else:
            continue

        if not filename or not download_url:
            continue

        # Filter by patterns if specified
        if file_patterns:
            if not any(p.lower() in filename.lower() for p in file_patterns):
                continue

        # Only download processed results, not raw files (too large)
        skip_extensions = [".raw", ".wiff", ".d", ".mzML", ".mzXML"]
        if any(filename.endswith(ext) for ext in skip_extensions):
            logger.info(f"Skipping raw file: {filename}")
            continue

        dest = dest_dir / filename
        try:
            download_file(download_url, dest)
            downloaded.append(dest)
        except Exception as e:
            logger.warning(f"Failed to download {filename}: {e}")

    return downloaded


def download_geo_supplementary(accession: str, file_patterns: list[str] | None = None) -> list[Path]:
    """Download supplementary files from GEO.

    Args:
        accession: GSE accession number.
        file_patterns: Optional filename filters.

    Returns:
        List of downloaded file paths.
    """
    dest_dir = DATA_DIR / accession
    dest_dir.mkdir(parents=True, exist_ok=True)

    # GEO FTP base
    prefix = accession[:7]  # GSE1234
    ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{accession}/suppl/"

    logger.info(f"Fetching GEO supplementary files: {accession}")

    try:
        response = requests.get(ftp_base, timeout=60)
        if response.status_code == 200:
            # Parse HTML directory listing for file links
            import re
            links = re.findall(r'href="([^"]+)"', response.text)
            data_links = [l for l in links if not l.startswith("?") and not l.startswith("/")]

            downloaded = []
            for link in data_links:
                if file_patterns and not any(p.lower() in link.lower() for p in file_patterns):
                    continue

                url = ftp_base + link
                dest = dest_dir / link
                try:
                    download_file(url, dest)
                    downloaded.append(dest)
                except Exception as e:
                    logger.warning(f"Failed to download {link}: {e}")

            return downloaded
    except Exception as e:
        logger.warning(f"Failed to access GEO FTP: {e}")

    return []


def download_stomics(dataset_id: str) -> list[Path]:
    """Download from STOmicsDB.

    Note: STOmicsDB may require manual download for some datasets.
    This function provides instructions and downloads what's programmatically accessible.
    """
    dest_dir = DATA_DIR / dataset_id
    dest_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"STOmicsDB dataset: {dataset_id}")
    logger.info(f"  Browse at: https://db.cngb.org/stomics/datasets/{dataset_id}")

    # STOmicsDB API
    api_url = f"https://db.cngb.org/stomics/api/datasets/{dataset_id}/"

    try:
        response = requests.get(api_url, timeout=60)
        if response.status_code == 200:
            data = response.json()
            # Save metadata
            meta_path = dest_dir / "metadata.json"
            with open(meta_path, "w") as f:
                json.dump(data, f, indent=2)
            logger.info(f"Saved metadata to {meta_path}")

            # Try to find download links in the response
            if "files" in data:
                downloaded = []
                for file_info in data["files"]:
                    url = file_info.get("url", file_info.get("download_url", ""))
                    name = file_info.get("name", file_info.get("filename", ""))
                    if url and name:
                        dest = dest_dir / name
                        try:
                            download_file(url, dest)
                            downloaded.append(dest)
                        except Exception as e:
                            logger.warning(f"Failed: {name}: {e}")
                return downloaded
    except Exception as e:
        logger.info(f"STOmicsDB API not accessible: {e}")

    # Create instructions file for manual download
    instructions = dest_dir / "DOWNLOAD_INSTRUCTIONS.md"
    instructions.write_text(f"""# Manual Download Required: {dataset_id}

## STOmicsDB Stereo-seq Aging Atlas

1. Visit: https://db.cngb.org/stomics/datasets/{dataset_id}
2. Download the GEM/GEF files for liver tissue (young vs old)
3. Place files in: {dest_dir}/

### Expected files:
- Liver_young.gem.gz (or similar)
- Liver_old.gem.gz (or similar)
- metadata.tsv

### Alternative: CNGB FTP
Try: https://ftp.cngb.org/pub/SciRAID/stomics/{dataset_id}/

### Alternative: Use the paper's GEO accession
The paper "Spatial transcriptomic landscape unveils immunoglobin-associated
senescence as a hallmark of aging" (Cell, 2024) may have GEO-deposited data.
Search GEO for the paper title.
""")

    logger.info(f"Manual download instructions written to {instructions}")
    return [instructions]


def download_elife_supplement(dest_dir: Path) -> list[Path]:
    """Download processed protein matrix from eLife PXD011967 supplement."""
    dest_dir.mkdir(parents=True, exist_ok=True)

    # The eLife paper (49874) has processed data as supplementary
    # Direct supplement download
    urls = [
        (
            "https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNDk4NzQvZWxpZmUtNDk4NzQtc3VwcDEtdjIueGxzeA--/elife-49874-supp1-v2.xlsx",
            "elife-49874-supp1-human_muscle_proteomics.xlsx",
        ),
    ]

    downloaded = []
    for url, name in urls:
        dest = dest_dir / name
        try:
            download_file(url, dest)
            downloaded.append(dest)
        except Exception as e:
            logger.warning(f"Failed to download eLife supplement: {e}")
            # Provide alternative
            alt_instructions = dest_dir / "DOWNLOAD_ELIFE.md"
            alt_instructions.write_text(f"""# Download PXD011967 Processed Data

1. Visit: https://elifesciences.org/articles/49874
2. Scroll to "Additional files" section
3. Download Supplementary File 1 (protein quantification matrix)
4. Place in: {dest_dir}/

Alternative: ProteomeXchange
- Visit: https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD011967
- Download processed result files
""")
            downloaded.append(alt_instructions)

    return downloaded


def main():
    """Download all datasets for the HBAM pipeline."""
    logger.remove()
    logger.add(sys.stderr, level="INFO", format="<green>{time:HH:mm:ss}</green> | {message}")

    print("=" * 60)
    print("HBAM Pipeline - Public Data Download")
    print("=" * 60)

    all_files = []

    # 1. PXD047296 - Ten Organs Proteome Atlas (liver + muscle proteomics)
    print("\n[1/3] PXD047296: Ten Organs Proteome Atlas")
    print("    Mouse liver + muscle proteomics, 4 timepoints")
    files = download_pride_files(
        "PXD047296",
        file_patterns=["proteingroups", "protein_groups", "report", "results", ".txt", ".tsv", ".csv", ".xlsx"],
    )
    all_files.extend(files)

    # 2. STDS0000247 - Stereo-seq Aging Atlas
    print("\n[2/3] STDS0000247: Stereo-seq Aging Atlas")
    print("    Mouse liver spatial transcriptomics, young vs old")
    files = download_stomics("STDS0000247")
    all_files.extend(files)

    # 3. PXD011967 - Human Muscle Aging Proteomics
    print("\n[3/3] PXD011967: Human Muscle Aging Proteomics")
    print("    58 subjects, 5 age strata, skeletal muscle")
    files = download_pride_files(
        "PXD011967",
        file_patterns=["proteingroups", "protein_groups", "report", "results", ".txt", ".tsv", ".csv", ".xlsx"],
    )
    # Also try eLife supplement
    elife_files = download_elife_supplement(DATA_DIR / "PXD011967")
    all_files.extend(files)
    all_files.extend(elife_files)

    # Summary
    print("\n" + "=" * 60)
    print(f"Download complete: {len(all_files)} files")
    for f in all_files:
        size = f.stat().st_size / 1e6 if f.exists() else 0
        print(f"  {f} ({size:.1f} MB)")

    print("\nNext steps:")
    print("  1. Check data/raw/ for downloaded files")
    print("  2. Update config/default.yaml with file paths")
    print("  3. Run: python -m hbam run --config config/default.yaml")


if __name__ == "__main__":
    main()
