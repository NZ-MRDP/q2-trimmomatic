"""Tests for local split/apply/combine helpers."""

import tempfile
import unittest
from pathlib import Path

import pandas as pd
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)

from q2_trimmomatic._pipelines import (
    collate_trimmed_paired,
    collate_trimmed_single,
    split_paired_samples,
    split_single_samples,
)


class PipelineHelperTests(unittest.TestCase):
    """Cover partitioning and collation helpers used by parallel pipelines."""

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

    def _make_trimmed_partition(self, filenames):
        """Create a minimal trimmed-output partition directory."""
        partition = CasavaOneEightSingleLanePerSampleDirFmt()
        for filename in filenames:
            Path(partition.path, filename).touch()
        return partition

    def test_split_single_samples_creates_local_partitions(self):
        """Single-end partitions should retain all samples and metadata."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sequences = self._make_single_reads_dir(tmpdir, ["a", "b", "c"])

            partitions = split_single_samples(sequences, num_partitions=2)

            self.assertEqual(set(partitions.keys()), {"0", "1"})

            observed_samples = set()
            for partition in partitions.values():
                self.assertTrue((partition.path / "metadata.yml").exists())
                partition_df = partition.manifest.view(pd.DataFrame)
                observed_samples.update(partition_df.index)

            self.assertEqual(observed_samples, {"a", "b", "c"})

    def test_split_paired_samples_caps_partition_count_to_sample_count(self):
        """Paired-end partitioning should not create empty artifacts."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paired_sequences = self._make_paired_reads_dir(tmpdir, ["a", "b"])

            partitions = split_paired_samples(paired_sequences, num_partitions=10)

            self.assertEqual(set(partitions.keys()), {"0", "1"})
            for partition in partitions.values():
                self.assertTrue((partition.path / "metadata.yml").exists())
                partition_df = partition.manifest.view(pd.DataFrame)
                self.assertEqual(len(partition_df.index), 1)

    def test_collate_trimmed_single_combines_all_partition_files(self):
        """Single-end collation should merge files from all partitions."""
        partitions = [
            self._make_trimmed_partition(["a_R1.fastq.gz"]),
            self._make_trimmed_partition(["b_R1.fastq.gz"]),
        ]

        result = collate_trimmed_single(partitions)

        self.assertTrue((result.path / "a_R1.fastq.gz").exists())
        self.assertTrue((result.path / "b_R1.fastq.gz").exists())

    def test_collate_trimmed_paired_accepts_collections(self):
        """Paired-end collation should accept mapping-like collections."""
        partitions = {
            "0": self._make_trimmed_partition(["a_R1.fastq.gz", "a_R2.fastq.gz"]),
            "1": self._make_trimmed_partition(["b_R1.fastq.gz", "b_R2.fastq.gz"]),
        }

        result = collate_trimmed_paired(partitions)

        self.assertTrue((result.path / "a_R1.fastq.gz").exists())
        self.assertTrue((result.path / "a_R2.fastq.gz").exists())
        self.assertTrue((result.path / "b_R1.fastq.gz").exists())
        self.assertTrue((result.path / "b_R2.fastq.gz").exists())
