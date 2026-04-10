"""CLI entry point for the HBAM pipeline."""

from __future__ import annotations

from pathlib import Path

import typer
from loguru import logger

app = typer.Typer(
    name="hbam",
    help="HBAM: High Burden Aging Muscle multi-omics pipeline",
    add_completion=False,
)


@app.command()
def run(
    config: Path = typer.Option(..., "--config", "-c", help="Path to YAML config file"),
    output_dir: Path = typer.Option(None, "--output-dir", "-o", help="Override output directory"),
    seed: int = typer.Option(None, "--seed", "-s", help="Override random seed"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate config without running"),
) -> None:
    """Run the full HBAM pipeline."""
    from hbam.config import load_config

    cfg = load_config(config)

    if seed is not None:
        cfg.seed = seed
    if output_dir is not None:
        cfg.data.output_dir = output_dir

    if dry_run:
        typer.echo(f"Config validated successfully:")
        typer.echo(f"  Seed: {cfg.seed}")
        typer.echo(f"  Modalities: {len(cfg.data.modalities)}")
        typer.echo(f"  Output: {cfg.data.output_dir}")
        typer.echo(f"  Embedding: {cfg.integration.embedding_method}")
        typer.echo(f"  Gene sets: {cfg.engine.gene_set_source}")
        return

    from hbam.pipeline import run_pipeline

    result = run_pipeline(cfg)

    typer.echo(f"\nPipeline complete!")
    typer.echo(f"  Figures: {len(result.figure_paths)}")
    typer.echo(f"  Tables: {len(result.table_paths)}")

    if result.scores is not None:
        typer.echo(f"  Samples scored: {len(result.scores)}")
        typer.echo(f"  HBAM range: [{result.scores['hbam_score'].min():.4f}, {result.scores['hbam_score'].max():.4f}]")


@app.command()
def qc(
    config: Path = typer.Option(..., "--config", "-c", help="Path to YAML config file"),
) -> None:
    """Run only data loading and QC."""
    from hbam.config import load_config
    from hbam.pipeline import _load_all_modalities, _process_modalities, set_seeds
    from hbam.utils.logging import setup_logging

    cfg = load_config(config)
    setup_logging(log_level=cfg.log_level)
    set_seeds(cfg.seed)

    modalities = _load_all_modalities(cfg)
    processed = _process_modalities(modalities, cfg)

    for name, adata in processed.items():
        typer.echo(f"{name}: {adata.n_obs} samples x {adata.n_vars} features")


@app.command()
def figures(
    config: Path = typer.Option(..., "--config", "-c", help="Path to YAML config file"),
) -> None:
    """Regenerate figures from saved intermediate results."""
    typer.echo("Figure regeneration requires intermediate results in results/intermediate/")
    typer.echo("Run the full pipeline first with: hbam run --config <config.yaml>")


@app.command()
def integrate(
    config: Path = typer.Option(..., "--config", "-c", help="Path to YAML config file"),
) -> None:
    """Run only integration (assumes processed data exists)."""
    from hbam.config import load_config
    from hbam.pipeline import _load_all_modalities, _process_modalities, _integrate_modalities, set_seeds
    from hbam.utils.logging import setup_logging

    cfg = load_config(config)
    setup_logging(log_level=cfg.log_level)
    set_seeds(cfg.seed)

    modalities = _load_all_modalities(cfg)
    if not modalities:
        from hbam.pipeline import _generate_demo_data
        modalities = _generate_demo_data(cfg.seed)

    processed = _process_modalities(modalities, cfg)
    mudata = _integrate_modalities(processed, cfg)

    for name, mod in mudata.mod.items():
        typer.echo(f"{name}: {mod.n_obs} samples x {mod.n_vars} genes")

    if "pca_factors" in mudata.obsm:
        typer.echo(f"Latent factors: {mudata.obsm['pca_factors'].shape[1]} components")


@app.command()
def score(
    config: Path = typer.Option(..., "--config", "-c", help="Path to YAML config file"),
) -> None:
    """Run only HBAM scoring (assumes integration done)."""
    from hbam.config import load_config
    from hbam.pipeline import (
        _load_all_modalities, _process_modalities, _integrate_modalities,
        _run_hbam_engine, _run_translation, set_seeds,
    )
    from hbam.utils.logging import setup_logging

    cfg = load_config(config)
    setup_logging(log_level=cfg.log_level)
    set_seeds(cfg.seed)

    modalities = _load_all_modalities(cfg)
    if not modalities:
        from hbam.pipeline import _generate_demo_data
        modalities = _generate_demo_data(cfg.seed)

    processed = _process_modalities(modalities, cfg)
    mudata = _integrate_modalities(processed, cfg)
    weights, gene_sets, enrichment = _run_hbam_engine(mudata, cfg)
    scores = _run_translation(mudata, weights, cfg)

    typer.echo(f"Scored {len(scores)} samples")
    typer.echo(f"HBAM range: [{scores['hbam_score'].min():.4f}, {scores['hbam_score'].max():.4f}]")


if __name__ == "__main__":
    app()
