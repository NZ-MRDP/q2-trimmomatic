"""QIIME 2 plugin for adapter trimming and quality control with Trimmomatic."""

from ._pipelines import trim_paired, trim_single

__version__ = "0.1.0"

__all__ = ["trim_paired", "trim_single"]
