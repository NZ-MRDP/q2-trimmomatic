"""QIIME 2 plugin for adapter trimming and quality control with Trimmomatic."""

from ._trimmomatic import trim_paired, trim_single
from ._pipelines import (
    collate_trimmed_paired,
    collate_trimmed_single,
    split_paired_samples,
    split_single_samples,
    trim_paired_parallel,
    trim_single_parallel,
)

__version__ = "0.1.0"

__all__ = [
    "trim_paired",
    "trim_single",
    "split_paired_samples",
    "split_single_samples",
    "collate_trimmed_paired",
    "collate_trimmed_single",
    "trim_paired_parallel",
    "trim_single_parallel",
]
