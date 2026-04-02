"""Tests for Trimmomatic command construction and internal parallel collation."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)

from q2_trimmomatic._trimmomatic import (
    _collate_partition,
    trim_paired,
    trim_single,
)


class TrimMethodTests(unittest.TestCase):
    """Cover CLI construction and internal result collation."""

    def _make_single_reads_dir(self, tmpdir, sample_ids):
        """Create a minimal single-end per-sample directory format."""
        manifest_rows = []
        for sample_id in sample_ids:
            filename = f"{sample_id}_S1_L001_R1_001.fastq.gz"
            Path(tmpdir, filename).touch()
            manifest_rows.append(f"{sample_id},{filename},forward")

        Path(tmpdir, "MANIFEST").write_text(
            "sample-id,filename,direction\n" + "\n".join(manifest_rows) + "\n"
        )
        Path(tmpdir, "metadata.yml").write_text("phred-offset: 33\n")

        return SingleLanePerSampleSingleEndFastqDirFmt(tmpdir, mode="r")

    def _make_paired_reads_dir(self, tmpdir, sample_ids):
        """Create a minimal paired-end per-sample directory format."""
        manifest_rows = []
        for sample_id in sample_ids:
            fwd = f"{sample_id}_S1_L001_R1_001.fastq.gz"
            rev = f"{sample_id}_S1_L001_R2_001.fastq.gz"
            Path(tmpdir, fwd).touch()
            Path(tmpdir, rev).touch()
            manifest_rows.append(f"{sample_id},{fwd},forward")
            manifest_rows.append(f"{sample_id},{rev},reverse")

        Path(tmpdir, "MANIFEST").write_text(
            "sample-id,filename,direction\n" + "\n".join(manifest_rows) + "\n"
        )
        Path(tmpdir, "metadata.yml").write_text("phred-offset: 33\n")

        return SingleLanePerSamplePairedEndFastqDirFmt(tmpdir, mode="r")

    def _subprocess_side_effect(self, cmd, check):
        """Create dummy output files in the paths Trimmomatic would write."""
        del check
        mode_idx = cmd.index("PE") if "PE" in cmd else cmd.index("SE")
        cursor = mode_idx + 1
        if cmd[cursor] == "-threads":
            cursor += 2

        if cmd[mode_idx] == "PE":
            output_paths = cmd[cursor + 2:cursor + 6]
        else:
            output_paths = [cmd[cursor + 1]]

        for output_path in output_paths:
            Path(output_path).touch()

    def test_collate_partition_copies_files(self):
        """Partition collation should copy all files into the destination."""
        source = CasavaOneEightSingleLanePerSampleDirFmt()
        destination = CasavaOneEightSingleLanePerSampleDirFmt()
        Path(source.path / "sample.fastq.gz").write_text("data")

        _collate_partition(source, destination)

        self.assertTrue((destination.path / "sample.fastq.gz").exists())

    @patch("q2_trimmomatic._trimmomatic.subprocess.run")
    def test_trim_single_omits_threads_flag_for_auto(self, mock_run):
        """Automatic thread selection should omit Trimmomatic's thread flag."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sequences = self._make_single_reads_dir(tmpdir, ["sampleA"])

            trim_single(sequences=sequences, threads="auto")

        cmd = mock_run.call_args[0][0]
        self.assertNotIn("-threads", cmd)

    @patch("q2_trimmomatic._trimmomatic.subprocess.run")
    def test_trim_paired_includes_threads_flag_for_explicit_value(self, mock_run):
        """Explicit thread counts should be forwarded to Trimmomatic."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paired_sequences = self._make_paired_reads_dir(tmpdir, ["sampleA"])

            trim_paired(paired_sequences=paired_sequences, threads=4)

        cmd = mock_run.call_args[0][0]
        self.assertIn("-threads", cmd)
        self.assertIn("4", cmd)

    @patch("q2_trimmomatic._trimmomatic.subprocess.run")
    def test_trim_single_collates_outputs_from_multiple_jobs(self, mock_run):
        """Single-end trimming should collate outputs from multiple jobs."""
        mock_run.side_effect = self._subprocess_side_effect
        with tempfile.TemporaryDirectory() as tmpdir:
            sequences = self._make_single_reads_dir(tmpdir, ["sampleA", "sampleB"])

            result = trim_single(sequences=sequences, n_jobs=2)

        self.assertEqual(mock_run.call_count, 2)
        self.assertTrue((result.path / "sampleA_S1_L001_R1_001.fastq.gz").exists())
        self.assertTrue((result.path / "sampleB_S1_L001_R1_001.fastq.gz").exists())

    @patch("q2_trimmomatic._trimmomatic.subprocess.run")
    def test_trim_paired_collates_outputs_from_multiple_jobs(self, mock_run):
        """Paired-end trimming should collate outputs from multiple jobs."""
        mock_run.side_effect = self._subprocess_side_effect
        with tempfile.TemporaryDirectory() as tmpdir:
            paired_sequences = self._make_paired_reads_dir(tmpdir, ["sampleA", "sampleB"])

            paired, unpaired_fwd, unpaired_rev = trim_paired(
                paired_sequences=paired_sequences,
                n_jobs=2,
            )

        self.assertEqual(mock_run.call_count, 2)
        self.assertTrue((paired.path / "sampleA_S1_L001_R1_001.fastq.gz").exists())
        self.assertTrue((paired.path / "sampleA_S1_L001_R2_001.fastq.gz").exists())
        self.assertTrue((paired.path / "sampleB_S1_L001_R1_001.fastq.gz").exists())
        self.assertTrue((paired.path / "sampleB_S1_L001_R2_001.fastq.gz").exists())
        self.assertTrue((unpaired_fwd.path / "sampleA_S1_L001_R1_001.fastq.gz").exists())
        self.assertTrue((unpaired_fwd.path / "sampleB_S1_L001_R1_001.fastq.gz").exists())
        self.assertTrue((unpaired_rev.path / "sampleA_S1_L001_R2_001.fastq.gz").exists())
        self.assertTrue((unpaired_rev.path / "sampleB_S1_L001_R2_001.fastq.gz").exists())
