"""Parallel trimmomatic split/apply/combine helpers and pipelines."""

import os
from pathlib import Path

import pandas as pd
from qiime2.util import duplicate

from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)


def split_paired_samples(
    paired_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
    num_partitions: int | None = None,
) -> SingleLanePerSamplePairedEndFastqDirFmt:
    """Split paired-end samples into a collection of smaller artifacts."""
    return _split_sequences(
        sequences=paired_sequences,
        output_factory=SingleLanePerSamplePairedEndFastqDirFmt,
        num_partitions=num_partitions,
        num_reads_per_sample=2,
    )


def split_single_samples(
    sequences: SingleLanePerSampleSingleEndFastqDirFmt,
    num_partitions: int | None = None,
) -> SingleLanePerSampleSingleEndFastqDirFmt:
    """Split single-end samples into a collection of smaller artifacts."""
    return _split_sequences(
        sequences=sequences,
        output_factory=SingleLanePerSampleSingleEndFastqDirFmt,
        num_partitions=num_partitions,
        num_reads_per_sample=1,
    )


def collate_trimmed_paired(
    trimmed: CasavaOneEightSingleLanePerSampleDirFmt,
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """Collate trimmed paired-end sequence partitions into a single artifact."""
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    for partition in _iter_collection_values(trimmed):
        for fp in partition.path.iterdir():
            duplicate(fp, result.path / fp.name)
    return result


def collate_trimmed_single(
    trimmed: CasavaOneEightSingleLanePerSampleDirFmt,
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """Collate trimmed single-end sequence partitions into a single artifact."""
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    for partition in _iter_collection_values(trimmed):
        for fp in partition.path.iterdir():
            duplicate(fp, result.path / fp.name)
    return result


def trim_paired_parallel(
    ctx,
    paired_sequences,
    num_partitions=None,
    threads=1,
    adapter_file="NexteraPE-PE.fa",
    leading=3,
    trailing=3,
    sliding_window_size=4,
    sliding_window_quality=15,
    min_length=100,
    head_crop=0,
    crop=0,
    seed_mismatches=2,
    palindrome_clip_threshold=30,
    simple_clip_threshold=10,
):
    """Parallel adapter trimming and quality control for paired-end reads."""
    _trim_paired = ctx.get_action("trimmomatic", "trim_paired")
    _split = ctx.get_action("trimmomatic", "split_paired_samples")
    _collate_pe = ctx.get_action("trimmomatic", "collate_trimmed_paired")
    _collate_se = ctx.get_action("trimmomatic", "collate_trimmed_single")

    (partitioned_seqs,) = _split(paired_sequences, num_partitions=num_partitions)

    pe_trimmed_parts = []
    unpaired_fwd_parts = []
    unpaired_rev_parts = []

    for seqs in partitioned_seqs.values():
        pe, uf, ur = _trim_paired(
            paired_sequences=seqs,
            threads=threads,
            adapter_file=adapter_file,
            leading=leading,
            trailing=trailing,
            sliding_window_size=sliding_window_size,
            sliding_window_quality=sliding_window_quality,
            min_length=min_length,
            head_crop=head_crop,
            crop=crop,
            seed_mismatches=seed_mismatches,
            palindrome_clip_threshold=palindrome_clip_threshold,
            simple_clip_threshold=simple_clip_threshold,
        )
        pe_trimmed_parts.append(pe)
        unpaired_fwd_parts.append(uf)
        unpaired_rev_parts.append(ur)

    (paired_end_trimmed,) = _collate_pe(pe_trimmed_parts)
    (unpaired_fwd,) = _collate_se(unpaired_fwd_parts)
    (unpaired_rev,) = _collate_se(unpaired_rev_parts)

    return paired_end_trimmed, unpaired_fwd, unpaired_rev


def trim_single_parallel(
    ctx,
    sequences,
    num_partitions=None,
    threads=1,
    adapter_file="TruSeq3-SE.fa",
    leading=3,
    trailing=3,
    sliding_window_size=4,
    sliding_window_quality=15,
    min_length=100,
    head_crop=0,
    crop=0,
    seed_mismatches=2,
    simple_clip_threshold=10,
):
    """Parallel adapter trimming and quality control for single-end reads."""
    _trim_single = ctx.get_action("trimmomatic", "trim_single")
    _split = ctx.get_action("trimmomatic", "split_single_samples")
    _collate_se = ctx.get_action("trimmomatic", "collate_trimmed_single")

    (partitioned_seqs,) = _split(sequences, num_partitions=num_partitions)

    trimmed_parts = []
    for seqs in partitioned_seqs.values():
        (trimmed,) = _trim_single(
            sequences=seqs,
            threads=threads,
            adapter_file=adapter_file,
            leading=leading,
            trailing=trailing,
            sliding_window_size=sliding_window_size,
            sliding_window_quality=sliding_window_quality,
            min_length=min_length,
            head_crop=head_crop,
            crop=crop,
            seed_mismatches=seed_mismatches,
            simple_clip_threshold=simple_clip_threshold,
        )
        trimmed_parts.append(trimmed)

    (collated_trimmed,) = _collate_se(trimmed_parts)
    return collated_trimmed


def _split_sequences(sequences, output_factory, num_partitions, num_reads_per_sample):
    """Create self-contained partition artifacts from an input sequence artifact."""
    df = sequences.manifest.view(pd.DataFrame)
    sample_ids = list(df.index)
    partition_ids = _partition_ids(sample_ids, num_partitions)
    result = {}

    for idx, sample_id_partition in enumerate(partition_ids):
        partition = output_factory()
        _copy_if_exists(sequences.path / "metadata.yml", partition.path / "metadata.yml")

        manifest_rows = []
        partition_df = df.loc[sample_id_partition]
        for sample_id, row in partition_df.iterrows():
            filenames = row.tolist()
            directions = _directions_for_sample(num_reads_per_sample)
            for filename, direction in zip(filenames, directions, strict=True):
                source = Path(filename)
                destination_name = os.path.basename(filename)
                duplicate(source, partition.path / destination_name)
                manifest_rows.append((sample_id, destination_name, direction))

        _write_manifest(partition.path / "MANIFEST", manifest_rows)
        result[str(idx)] = partition

    return result


def _partition_ids(sample_ids, num_partitions):
    """Spread sample ids across a bounded number of non-empty partitions."""
    if not sample_ids:
        return []

    if num_partitions is None:
        num_partitions = len(sample_ids)

    num_partitions = max(1, min(num_partitions, len(sample_ids)))
    return [sample_ids[idx::num_partitions] for idx in range(num_partitions)]


def _directions_for_sample(num_reads_per_sample):
    """Return manifest direction labels for one or two reads per sample."""
    if num_reads_per_sample == 1:
        return ("forward",)

    return ("forward", "reverse")


def _write_manifest(path, rows):
    """Write a per-sample FASTQ manifest for a partition artifact."""
    with open(path, "w") as fh:
        fh.write("sample-id,filename,direction\n")
        for sample_id, filename, direction in rows:
            fh.write(f"{sample_id},{filename},{direction}\n")


def _copy_if_exists(source, destination):
    """Duplicate a file into the destination path when it exists."""
    if source.exists():
        duplicate(source, destination)


def _iter_collection_values(collection):
    """Iterate over list-like or mapping-like QIIME collections."""
    if hasattr(collection, "values"):
        return collection.values()

    return collection
