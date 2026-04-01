"""Parallel trimmomatic pipelines using the QIIME 2 pipeline/ctx pattern."""

from qiime2.util import duplicate

from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
)


def collate_trimmed_paired(
    trimmed: CasavaOneEightSingleLanePerSampleDirFmt,
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """Collate trimmed paired-end sequence partitions into a single artifact."""
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    for partition in trimmed:
        for fp in partition.path.iterdir():
            duplicate(fp, result.path / fp.name)
    return result


def collate_trimmed_single(
    trimmed: CasavaOneEightSingleLanePerSampleDirFmt,
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """Collate trimmed single-end sequence partitions into a single artifact."""
    result = CasavaOneEightSingleLanePerSampleDirFmt()
    for partition in trimmed:
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
    """Parallel adapter trimming and quality control for paired-end reads.

    Partitions samples across workers using ``demux.partition_samples_paired``,
    trims each partition independently with ``trimmomatic.trim_paired``, then
    collates the results back into single artifacts.  When a parallel execution
    backend (e.g. parsl) is configured, each partition runs concurrently.

    Parameters
    ----------
    ctx : qiime2.sdk.Context
        QIIME 2 pipeline context used to look up registered actions.
    paired_sequences : SampleData[PairedEndSequencesWithQuality]
        Paired-end FASTQ sequence data to trim.
    num_partitions : int, optional
        Number of partitions to split samples into. Defaults to one partition
        per sample.
    adapter_file : str, optional
        Bundled adapter FASTA file for ILLUMINACLIP. Default is
        "NexteraPE-PE.fa".
    leading : int, optional
        Minimum quality to keep a base at the 5' end (LEADING). Default 3.
    trailing : int, optional
        Minimum quality to keep a base at the 3' end (TRAILING). Default 3.
    sliding_window_size : int, optional
        Sliding window width for quality averaging (SLIDINGWINDOW). Default 4.
    sliding_window_quality : int, optional
        Minimum average quality within the sliding window. Default 15.
    min_length : int, optional
        Minimum read length after trimming (MINLEN). Default 100.
    head_crop : int, optional
        Bases to remove from the 5' end (HEADCROP). Set to 0 to disable.
    crop : int, optional
        Cut reads to this length from the 3' end (CROP). Set to 0 to disable.
    seed_mismatches : int, optional
        Maximum mismatches in the ILLUMINACLIP seed. Default 2.
    palindrome_clip_threshold : int, optional
        Threshold for palindrome adapter clipping (PE only). Default 30.
    simple_clip_threshold : int, optional
        Threshold for simple adapter clipping. Default 10.

    Returns
    -------
    SampleData[PairedEndSequencesWithQuality]
        Trimmed paired-end reads where both mates passed filtering.
    SampleData[SequencesWithQuality]
        Trimmed forward reads whose reverse mate failed filtering.
    SampleData[SequencesWithQuality]
        Trimmed reverse reads whose forward mate failed filtering.

    """
    _trim_paired = ctx.get_action("trimmomatic", "trim_paired")
    _partition = ctx.get_action("demux", "partition_samples_paired")
    _collate_pe = ctx.get_action("trimmomatic", "collate_trimmed_paired")
    _collate_se = ctx.get_action("trimmomatic", "collate_trimmed_single")

    (partitioned_seqs,) = _partition(paired_sequences, num_partitions)

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
    """Parallel adapter trimming and quality control for single-end reads.

    Partitions samples across workers using ``demux.partition_samples_single``,
    trims each partition independently with ``trimmomatic.trim_single``, then
    collates the results back into a single artifact.  When a parallel
    execution backend (e.g. parsl) is configured, each partition runs
    concurrently.

    Parameters
    ----------
    ctx : qiime2.sdk.Context
        QIIME 2 pipeline context used to look up registered actions.
    sequences : SampleData[SequencesWithQuality]
        Single-end FASTQ sequence data to trim.
    num_partitions : int, optional
        Number of partitions to split samples into. Defaults to one partition
        per sample.
    adapter_file : str, optional
        Bundled adapter FASTA file for ILLUMINACLIP. Default is
        "TruSeq3-SE.fa".
    leading : int, optional
        Minimum quality to keep a base at the 5' end (LEADING). Default 3.
    trailing : int, optional
        Minimum quality to keep a base at the 3' end (TRAILING). Default 3.
    sliding_window_size : int, optional
        Sliding window width for quality averaging (SLIDINGWINDOW). Default 4.
    sliding_window_quality : int, optional
        Minimum average quality within the sliding window. Default 15.
    min_length : int, optional
        Minimum read length after trimming (MINLEN). Default 100.
    head_crop : int, optional
        Bases to remove from the 5' end (HEADCROP). Set to 0 to disable.
    crop : int, optional
        Cut reads to this length from the 3' end (CROP). Set to 0 to disable.
    seed_mismatches : int, optional
        Maximum mismatches in the ILLUMINACLIP seed. Default 2.
    simple_clip_threshold : int, optional
        Threshold for simple adapter clipping. Default 10.

    Returns
    -------
    SampleData[SequencesWithQuality]
        Trimmed single-end reads that passed filtering.

    """
    _trim_single = ctx.get_action("trimmomatic", "trim_single")
    _partition = ctx.get_action("demux", "partition_samples_single")
    _collate_se = ctx.get_action("trimmomatic", "collate_trimmed_single")

    (partitioned_seqs,) = _partition(sequences, num_partitions)

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
