"""Structured logging with loguru for the HBAM pipeline."""

from __future__ import annotations

import sys
import time
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Any, Generator

from loguru import logger


def setup_logging(
    log_level: str = "INFO",
    log_file: Path | None = None,
    log_dir: Path = Path("results/logs"),
) -> None:
    """Configure loguru for the pipeline.

    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR).
        log_file: Explicit log file path. If None, auto-generated in log_dir.
        log_dir: Directory for auto-generated log files.
    """
    # Remove default handler
    logger.remove()

    # Stderr handler (human-readable)
    logger.add(
        sys.stderr,
        level=log_level,
        format="<green>{time:HH:mm:ss}</green> | <level>{level:<7}</level> | <cyan>{module}</cyan>:<cyan>{function}</cyan> | {message}",
        colorize=True,
    )

    # File handler (structured, for parsing)
    if log_file is None:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = log_dir / f"pipeline_{timestamp}.log"
    else:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)

    logger.add(
        str(log_file),
        level="DEBUG",
        format="{time:YYYY-MM-DD HH:mm:ss.SSS} | {level:<7} | {module}:{function}:{line} | {message}",
        rotation="50 MB",
        retention="7 days",
    )

    logger.info(f"Logging initialized | level={log_level} | file={log_file}")


@contextmanager
def log_step(name: str, **extra: Any) -> Generator[None, None, None]:
    """Context manager that logs entry/exit with elapsed time.

    Args:
        name: Name of the processing step.
        **extra: Additional key-value pairs to log.

    Usage:
        with log_step("normalization", modality="proteomics"):
            normalize(data)
    """
    extra_str = " | ".join(f"{k}={v}" for k, v in extra.items())
    prefix = f"STEP | {name}"
    if extra_str:
        prefix = f"{prefix} | {extra_str}"

    logger.info(f"{prefix} | status=STARTED")
    start = time.perf_counter()

    try:
        yield
        elapsed = time.perf_counter() - start
        logger.info(f"{prefix} | status=COMPLETED | elapsed={elapsed:.2f}s")
    except Exception as e:
        elapsed = time.perf_counter() - start
        logger.error(f"{prefix} | status=FAILED | elapsed={elapsed:.2f}s | error={type(e).__name__}: {e}")
        raise


def log_filter(
    action: str,
    before: int,
    after: int,
    reason: str,
    **extra: Any,
) -> None:
    """Log a data filtering event with structured metrics.

    Args:
        action: Name of the filter operation.
        before: Number of items before filtering.
        after: Number of items after filtering.
        reason: Why items were removed.
        **extra: Additional context.
    """
    removed = before - after
    pct = (removed / before * 100) if before > 0 else 0.0
    extra_str = " | ".join(f"{k}={v}" for k, v in extra.items())
    msg = f"FILTER | action={action} | before={before} | after={after} | removed={removed} ({pct:.1f}%) | reason={reason}"
    if extra_str:
        msg = f"{msg} | {extra_str}"
    logger.info(msg)


def log_metric(name: str, value: Any, **extra: Any) -> None:
    """Log a computed metric.

    Args:
        name: Metric name.
        value: Metric value.
        **extra: Additional context.
    """
    extra_str = " | ".join(f"{k}={v}" for k, v in extra.items())
    msg = f"METRIC | {name}={value}"
    if extra_str:
        msg = f"{msg} | {extra_str}"
    logger.info(msg)
