"""Tests for Trimmomatic command construction."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)

from q2_trimmomatic._trimmomatic import trim_paired, trim_single


class ThreadsHandlingTests(unittest.TestCase):
    """Cover translation of QIIME thread settings into CLI arguments."""

    def _make_single_reads_dir(self, tmpdir):
        """Create a minimal single-end per-sample directory format."""
        filename = "sampleA_S1_L001_R1_001.fastq.gz"
        Path(tmpdir, filename).touch()
        Path(tmpdir, "MANIFEST").write_text(
            "sample-id,filename,direction\n"
            f"sampleA,{filename},forward\n"
        )
        Path(tmpdir, "metadata.yml").write_text("phred-offset: 33\n")

        return SingleLanePerSampleSingleEndFastqDirFmt(tmpdir, mode="r")

    def _make_paired_reads_dir(self, tmpdir):
        """Create a minimal paired-end per-sample directory format."""
        fwd = "sampleA_S1_L001_R1_001.fastq.gz"
        rev = "sampleA_S1_L001_R2_001.fastq.gz"
        Path(tmpdir, fwd).touch()
        Path(tmpdir, rev).touch()
        Path(tmpdir, "MANIFEST").write_text(
            "sample-id,filename,direction\n"
            f"sampleA,{fwd},forward\n"
            f"sampleA,{rev},reverse\n"
        )
        Path(tmpdir, "metadata.yml").write_text("phred-offset: 33\n")

        return SingleLanePerSamplePairedEndFastqDirFmt(tmpdir, mode="r")

    @patch("q2_trimmomatic._trimmomatic.subprocess.run")
    def test_trim_single_omits_threads_flag_for_auto(self, mock_run):
        """Automatic thread selection should omit Trimmomatic's thread flag."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sequences = self._make_single_reads_dir(tmpdir)

            trim_single(sequences=sequences, threads="auto")

        cmd = mock_run.call_args[0][0]
        self.assertNotIn("-threads", cmd)

    @patch("q2_trimmomatic._trimmomatic.subprocess.run")
    def test_trim_paired_includes_threads_flag_for_explicit_value(self, mock_run):
        """Explicit thread counts should be forwarded to Trimmomatic."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paired_sequences = self._make_paired_reads_dir(tmpdir)

            trim_paired(paired_sequences=paired_sequences, threads=4)

        cmd = mock_run.call_args[0][0]
        self.assertIn("-threads", cmd)
        self.assertIn("4", cmd)
