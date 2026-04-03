"""Public QIIME 2 pipelines and internal split/combine helper actions."""

import os
from pathlib import Path

import pandas as pd
from qiime2.util import duplicate

from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)


def trim_paired(
    ctx,
    paired_sequences,
    num_partitions=None,
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
    threads=1,
):
    """Trim paired-end reads using an internal split/apply/combine pipeline."""
    split_action = ctx.get_action("trimmomatic", "_split_paired_samples")
    worker_action = ctx.get_action("trimmomatic", "_trim_paired_worker")
    combine_paired_action = ctx.get_action("trimmomatic", "_combine_paired_sequences")
    combine_single_action = ctx.get_action("trimmomatic", "_combine_single_sequences")

    pieces, = split_action(paired_sequences, num_partitions=num_partitions)

    paired_results = {}
    unpaired_fwd_results = {}
    unpaired_rev_results = {}
    for key, piece in pieces.items():
        paired_piece, unpaired_fwd_piece, unpaired_rev_piece = worker_action(
            piece,
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
            threads=threads,
        )
        paired_results[key] = paired_piece
        unpaired_fwd_results[key] = unpaired_fwd_piece
        unpaired_rev_results[key] = unpaired_rev_piece

    paired_end_trimmed, = combine_paired_action(paired_results)
    unpaired_fwd, = combine_single_action(unpaired_fwd_results)
    unpaired_rev, = combine_single_action(unpaired_rev_results)
    return paired_end_trimmed, unpaired_fwd, unpaired_rev


def trim_single(
    ctx,
    sequences,
    num_partitions=None,
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
    threads=1,
):
    """Trim single-end reads using an internal split/apply/combine pipeline."""
    split_action = ctx.get_action("trimmomatic", "_split_single_samples")
    worker_action = ctx.get_action("trimmomatic", "_trim_single_worker")
    combine_action = ctx.get_action("trimmomatic", "_combine_single_sequences")

    pieces, = split_action(sequences, num_partitions=num_partitions)

    partial_results = {}
    for key, piece in pieces.items():
        partial_results[key], = worker_action(
            piece,
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
            threads=threads,
        )

    final_result, = combine_action(partial_results)
    return final_result


def _split_paired_samples(
    paired_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
    num_partitions: int | None = None,
) -> dict[str, SingleLanePerSamplePairedEndFastqDirFmt]:
    """Split paired-end samples into internal partition artifacts."""
    return _split_sequences(
        sequences=paired_sequences,
        output_factory=SingleLanePerSamplePairedEndFastqDirFmt,
        num_partitions=num_partitions,
        num_reads_per_sample=2,
    )


def _split_single_samples(
    sequences: SingleLanePerSampleSingleEndFastqDirFmt,
    num_partitions: int | None = None,
) -> dict[str, SingleLanePerSampleSingleEndFastqDirFmt]:
    """Split single-end samples into internal partition artifacts."""
    return _split_sequences(
        sequences=sequences,
        output_factory=SingleLanePerSampleSingleEndFastqDirFmt,
        num_partitions=num_partitions,
        num_reads_per_sample=1,
    )


def _combine_paired_sequences(
    trimmed: dict[str, CasavaOneEightSingleLanePerSampleDirFmt],
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """Combine paired-end partition outputs into one artifact."""
    return _combine_sequences(trimmed)


def _combine_single_sequences(
    trimmed: dict[str, CasavaOneEightSingleLanePerSampleDirFmt],
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """Combine single-end partition outputs into one artifact."""
    return _combine_sequences(trimmed)


def _split_sequences(sequences, output_factory, num_partitions, num_reads_per_sample):
    """Create self-contained partition artifacts from a sequence artifact."""
    df = sequences.manifest.view(pd.DataFrame)
    sample_ids = list(df.index)
    partitions = _partition_ids(sample_ids, num_partitions)
    result = {}

    for idx, sample_id_partition in enumerate(partitions):
        partition = output_factory()
        _copy_if_exists(sequences.path / "metadata.yml", partition.path / "metadata.yml")

        manifest_rows = []
        partition_df = df.loc[sample_id_partition]
        for sample_id, row in partition_df.iterrows():
            filenames = row.tolist()
            directions = _directions_for_sample(num_reads_per_sample)
            for filename, direction in zip(filenames, directions, strict=True):
                source = _resolve_manifest_path(sequences.path, filename)
                destination_name = os.path.basename(filename)
                duplicate(source, partition.path / destination_name)
                manifest_rows.append((sample_id, destination_name, direction))

        _write_manifest(partition.path / "MANIFEST", manifest_rows)
        result[str(idx)] = partition

    return result


def _combine_sequences(trimmed):
    """Collate partition files into a single output artifact."""
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    for partition in _iter_collection_values(trimmed):
        for fp in partition.path.iterdir():
            duplicate(fp, result.path / fp.name)
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


def _resolve_manifest_path(base_dir, filepath):
    """Resolve relative manifest entries against the source artifact path."""
    path = Path(filepath)
    if path.is_absolute():
        return path

    return base_dir / path


def _iter_collection_values(collection):
    """Iterate over mapping-like or sequence-like QIIME collections."""
    if hasattr(collection, "values"):
        return collection.values()

    return collection
