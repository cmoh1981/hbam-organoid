"""Pipeline configuration with Pydantic validation."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel, Field, field_validator


class QCConfig(BaseModel):
    """Quality control parameters."""
    missingness_threshold: float = Field(0.3, ge=0.0, le=1.0, description="Max fraction of missing values per feature")
    min_presence: float = Field(0.7, ge=0.0, le=1.0, description="Min fraction of samples a feature must be present in")
    variance_percentile: float = Field(10.0, ge=0.0, le=100.0, description="Remove features below this variance percentile")
    min_genes_per_sample: int = Field(200, ge=0, description="Min genes detected per sample (spatial)")
    min_samples_per_gene: int = Field(10, ge=0, description="Min samples a gene must appear in")
    max_mito_fraction: float = Field(0.2, ge=0.0, le=1.0, description="Max mitochondrial gene fraction (spatial)")


class NormConfig(BaseModel):
    """Normalization parameters."""
    proteomics_method: str = Field("median", description="Normalization method for proteomics: median, quantile")
    spatial_method: str = Field("scran", description="Normalization method for spatial: sctransform, scran, total")
    pseudocount: float = Field(1.0, ge=0.0, description="Pseudocount for log2 transform")
    log_transform: bool = Field(True, description="Apply log2 transform after normalization")


class ImputeConfig(BaseModel):
    """Missing value imputation parameters."""
    method: str = Field("auto", description="Imputation method: auto, mnar, mcar, none")
    mnar_method: str = Field("minprob", description="MNAR imputation: minprob, min")
    mnar_quantile: float = Field(0.01, ge=0.0, le=1.0, description="Quantile for MinProb imputation")
    mcar_method: str = Field("knn", description="MCAR imputation: knn, mean, median")
    knn_k: int = Field(5, ge=1, description="Number of neighbors for KNN imputation")


class IntegrationConfig(BaseModel):
    """Multi-omics integration parameters."""
    min_gene_overlap: int = Field(1000, ge=0, description="Minimum gene overlap required across modalities")
    scaling_method: str = Field("zscore", description="Feature scaling: zscore, minmax, none")
    embedding_method: str = Field("auto", description="Embedding method: auto, mofa, pca")
    n_factors: int = Field(15, ge=2, description="Number of latent factors/components")
    convergence_tolerance: float = Field(1e-4, gt=0, description="MOFA convergence tolerance")


class EngineConfig(BaseModel):
    """HBAM scoring engine parameters."""
    gene_set_source: str = Field("hybrid", description="Gene set source: literature, data_driven, hybrid, custom")
    custom_maturation_path: Optional[Path] = Field(None, description="Path to custom maturation gene set")
    custom_dysfunction_path: Optional[Path] = Field(None, description="Path to custom dysfunction gene set")
    refinement_method: str = Field("trajectory", description="Data-driven refinement: trajectory, clustering")
    bootstrap_iterations: int = Field(5, ge=1, description="Number of bootstrap iterations")
    bootstrap_fraction: float = Field(0.8, ge=0.1, le=1.0, description="Fraction of samples per bootstrap")
    bootstrap_min_spearman: float = Field(0.8, ge=0.0, le=1.0, description="Minimum Spearman correlation for validation")
    auto_flip_direction: bool = Field(True, description="Auto-flip HBAM sign if negatively correlated with age")


class TranslationConfig(BaseModel):
    """Translation layer parameters."""
    top_biomarkers: int = Field(20, ge=1, description="Number of top biomarkers to select")
    panel_sizes: list[int] = Field(default_factory=lambda: [5, 10, 15, 20, 30, 50], description="Panel sizes for optimization")
    effect_size_method: str = Field("cohens_d", description="Effect size: cohens_d, hedges_g")


class OutputConfig(BaseModel):
    """Output and figure parameters."""
    figure_formats: list[str] = Field(default_factory=lambda: ["pdf", "png"], description="Output figure formats")
    figure_dpi: int = Field(300, ge=72, description="DPI for raster output")
    single_col_width: float = Field(3.5, gt=0, description="Single column width in inches")
    double_col_width: float = Field(7.0, gt=0, description="Double column width in inches")
    font_family: str = Field("Arial", description="Font family for figures")
    font_size_label: int = Field(8, ge=4, description="Font size for axis labels")
    font_size_title: int = Field(10, ge=4, description="Font size for titles")
    colorblind_safe: bool = Field(True, description="Use colorblind-safe palette")
    top_n_heatmap: int = Field(50, ge=1, description="Top N genes for weight heatmap")
    top_n_volcano_labels: int = Field(10, ge=0, description="Top N genes to label on volcano plot")
    top_n_pathways: int = Field(20, ge=1, description="Top N pathways for enrichment plot")


class ModalityInput(BaseModel):
    """Configuration for a single data modality."""
    name: str = Field(..., description="Modality name (e.g., 'liver_proteomics')")
    path: Path = Field(..., description="Path to input data file")
    format: str = Field(..., description="File format: maxquant, diann, csv, tsv, gem")
    species: str = Field("mouse", description="Source species: mouse, human")
    modality_type: str = Field(..., description="Type: temporal, spatial, functional")
    time_col: Optional[str] = Field(None, description="Column name for timepoint (temporal only)")
    condition_col: Optional[str] = Field(None, description="Column name for condition/group")
    batch_col: Optional[str] = Field(None, description="Column name for batch (if batch correction needed)")


class DataConfig(BaseModel):
    """Data layer configuration."""
    modalities: list[ModalityInput] = Field(default_factory=list, description="List of input modality configurations")
    output_dir: Path = Field(Path("results"), description="Output directory")
    intermediate_dir: Path = Field(Path("results/intermediate"), description="Intermediate results directory")


class PipelineConfig(BaseModel):
    """Top-level pipeline configuration."""
    seed: int = Field(42, ge=0, description="Random seed for reproducibility")
    data: DataConfig = Field(default_factory=DataConfig)
    qc: QCConfig = Field(default_factory=QCConfig)
    normalization: NormConfig = Field(default_factory=NormConfig)
    imputation: ImputeConfig = Field(default_factory=ImputeConfig)
    integration: IntegrationConfig = Field(default_factory=IntegrationConfig)
    engine: EngineConfig = Field(default_factory=EngineConfig)
    translation: TranslationConfig = Field(default_factory=TranslationConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)
    log_level: str = Field("INFO", description="Logging level: DEBUG, INFO, WARNING, ERROR")
    log_file: Optional[Path] = Field(None, description="Log file path (auto-generated if None)")

    @field_validator("log_level")
    @classmethod
    def validate_log_level(cls, v: str) -> str:
        valid = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
        if v.upper() not in valid:
            raise ValueError(f"log_level must be one of {valid}")
        return v.upper()


def load_config(path: Path) -> PipelineConfig:
    """Load pipeline configuration from YAML file.

    Args:
        path: Path to YAML configuration file.

    Returns:
        Validated PipelineConfig object.

    Raises:
        FileNotFoundError: If config file doesn't exist.
        ValidationError: If config values are invalid.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with open(path) as f:
        raw = yaml.safe_load(f) or {}

    return PipelineConfig(**raw)


def save_config(config: PipelineConfig, path: Path) -> None:
    """Save pipeline configuration to YAML file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    data = config.model_dump(mode="json")
    # Convert Path objects to strings for YAML serialization
    _convert_paths(data)

    with open(path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def _convert_paths(d: dict) -> None:
    """Recursively convert Path objects to strings in a dict."""
    for k, v in d.items():
        if isinstance(v, dict):
            _convert_paths(v)
        elif isinstance(v, Path):
            d[k] = str(v)
