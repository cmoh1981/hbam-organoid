"""Validation gates for pipeline stage transitions."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from loguru import logger


class PipelineValidationError(Exception):
    """Raised when a validation gate fails."""

    def __init__(self, failures: list[ValidationResult]) -> None:
        self.failures = failures
        messages = [f"  - {f.name}: {f.message}" for f in failures]
        super().__init__(
            f"Validation gate failed with {len(failures)} error(s):\n" + "\n".join(messages)
        )


@dataclass
class ValidationResult:
    """Result of a single validation check."""
    name: str
    passed: bool
    message: str
    metrics: dict[str, Any] = field(default_factory=dict)


class ValidationGate:
    """Collects validation results and enforces pass/fail gates.

    Usage:
        gate = ValidationGate("post_qc")
        gate.add(validate_missingness(adata, 0.3))
        gate.add(validate_normalization(adata))
        gate.enforce()  # raises PipelineValidationError if any failed
    """

    def __init__(self, name: str) -> None:
        self.name = name
        self.results: list[ValidationResult] = []

    def add(self, result: ValidationResult) -> None:
        """Add a validation result."""
        self.results.append(result)
        status = "PASS" if result.passed else "FAIL"
        logger.info(f"VALIDATE | gate={self.name} | check={result.name} | {status} | {result.message}")

    def enforce(self) -> None:
        """Raise PipelineValidationError if any check failed."""
        failures = [r for r in self.results if not r.passed]
        if failures:
            raise PipelineValidationError(failures)
        logger.info(f"VALIDATE | gate={self.name} | ALL PASSED ({len(self.results)} checks)")

    @property
    def passed(self) -> bool:
        """Whether all checks passed."""
        return all(r.passed for r in self.results)

    @property
    def summary(self) -> dict[str, Any]:
        """Summary of all results."""
        return {
            "gate": self.name,
            "total": len(self.results),
            "passed": sum(1 for r in self.results if r.passed),
            "failed": sum(1 for r in self.results if not r.passed),
            "results": [
                {"name": r.name, "passed": r.passed, "message": r.message, "metrics": r.metrics}
                for r in self.results
            ],
        }


def validate_missingness(
    X: np.ndarray,
    threshold: float,
    name: str = "missingness",
) -> ValidationResult:
    """Check that fraction of missing values is below threshold.

    Args:
        X: Data matrix (samples x features).
        threshold: Maximum allowed fraction of missing values.
        name: Name for this check.

    Returns:
        ValidationResult with pass/fail and metrics.
    """
    if X.size == 0:
        return ValidationResult(
            name=name,
            passed=False,
            message="Empty matrix",
            metrics={"total_elements": 0},
        )

    nan_count = np.isnan(X).sum() if np.issubdtype(X.dtype, np.floating) else 0
    total = X.size
    fraction = nan_count / total
    passed = fraction <= threshold

    return ValidationResult(
        name=name,
        passed=passed,
        message=f"Missing fraction: {fraction:.4f} (threshold: {threshold})",
        metrics={"nan_count": int(nan_count), "total": total, "fraction": float(fraction)},
    )


def validate_dimensions(
    shape: tuple[int, ...],
    expected_shape: tuple[int | None, ...],
    name: str = "dimensions",
) -> ValidationResult:
    """Check matrix dimensions match expectations.

    Args:
        shape: Actual shape.
        expected_shape: Expected shape (None for any dimension).
        name: Name for this check.
    """
    mismatches = []
    for i, (actual, expected) in enumerate(zip(shape, expected_shape)):
        if expected is not None and actual != expected:
            mismatches.append(f"dim{i}: expected {expected}, got {actual}")

    passed = len(mismatches) == 0
    message = "Dimensions match" if passed else f"Dimension mismatch: {'; '.join(mismatches)}"

    return ValidationResult(
        name=name,
        passed=passed,
        message=message,
        metrics={"actual_shape": list(shape), "expected_shape": list(expected_shape)},
    )


def validate_gene_overlap(
    genes_a: set[str],
    genes_b: set[str],
    min_overlap: int = 1000,
    name: str = "gene_overlap",
) -> ValidationResult:
    """Check that gene overlap between two sets meets minimum.

    Args:
        genes_a: First gene set.
        genes_b: Second gene set.
        min_overlap: Minimum required overlap.
        name: Name for this check.
    """
    overlap = genes_a & genes_b
    n_overlap = len(overlap)
    smaller = min(len(genes_a), len(genes_b))
    fraction = n_overlap / smaller if smaller > 0 else 0.0
    passed = n_overlap >= min_overlap

    return ValidationResult(
        name=name,
        passed=passed,
        message=f"Overlap: {n_overlap} genes ({fraction:.1%} of smaller set, min required: {min_overlap})",
        metrics={
            "overlap_count": n_overlap,
            "set_a_size": len(genes_a),
            "set_b_size": len(genes_b),
            "fraction_of_smaller": float(fraction),
        },
    )


def validate_normalization(
    medians: np.ndarray,
    fold_range: float = 1.5,
    name: str = "normalization",
) -> ValidationResult:
    """Check that per-sample medians are within acceptable range.

    Args:
        medians: Array of per-sample median values.
        fold_range: Maximum fold-range between min and max median.
        name: Name for this check.
    """
    if len(medians) == 0:
        return ValidationResult(name=name, passed=False, message="No medians provided")

    med_min = float(np.min(medians))
    med_max = float(np.max(medians))

    if med_min <= 0:
        return ValidationResult(
            name=name,
            passed=False,
            message=f"Non-positive median detected: {med_min}",
            metrics={"min_median": med_min, "max_median": med_max},
        )

    actual_fold = med_max / med_min
    passed = actual_fold <= fold_range

    return ValidationResult(
        name=name,
        passed=passed,
        message=f"Median fold-range: {actual_fold:.2f} (max allowed: {fold_range})",
        metrics={"min_median": med_min, "max_median": med_max, "fold_range": float(actual_fold)},
    )
