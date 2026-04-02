"""Trimmomatic adapter trimming and quality control functions."""

from concurrent.futures import ThreadPoolExecutor
from importlib.resources import files
import os
import shutil
import subprocess

import pandas as pd
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)


def trim_paired(
    paired_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
    adapter_file: str = "NexteraPE-PE.fa",
    n_jobs: int = 1,
    leading: int = 3,
    trailing: int = 3,
    sliding_window_size: int = 4,
    sliding_window_quality: int = 15,
    min_length: int = 100,
    head_crop: int = 0,
    crop: int = 0,
    seed_mismatches: int = 2,
    palindrome_clip_threshold: int = 30,
    simple_clip_threshold: int = 10,
    threads: int | str = 1,
) -> (CasavaOneEightSingleLanePerSampleDirFmt,
      CasavaOneEightSingleLanePerSampleDirFmt,
      CasavaOneEightSingleLanePerSampleDirFmt):
    """Trim paired-end reads using Trimmomatic.

    Parameters
    ----------
    paired_sequences : SingleLanePerSamplePairedEndFastqDirFmt
        Paired-end FASTQ sequence data to trim.
    adapter_file : str, optional
        Name of the bundled adapter FASTA file to use for ILLUMINACLIP.
        Default is "NexteraPE-PE.fa".
    n_jobs : int, optional
        Number of samples to process concurrently. Each concurrent job
        launches a separate Trimmomatic process. Default is 1.
    leading : int, optional
        Minimum quality required to keep a base at the start of a read (LEADING).
        Default is 3.
    trailing : int, optional
        Minimum quality required to keep a base at the end of a read (TRAILING).
        Default is 3.
    sliding_window_size : int, optional
        Number of bases in the sliding window for quality averaging (SLIDINGWINDOW).
        Default is 4.
    sliding_window_quality : int, optional
        Minimum average quality required in the sliding window (SLIDINGWINDOW).
        Default is 15.
    min_length : int, optional
        Reads shorter than this length after trimming will be discarded (MINLEN).
        Default is 100.
    head_crop : int, optional
        Number of bases to remove from the start of each read (HEADCROP).
        Applied before adapter clipping. Set to 0 to disable. Default is 0.
    crop : int, optional
        Cut reads to this length by removing bases from the end (CROP).
        Applied after quality trimming but before MINLEN. Set to 0 to disable.
        Default is 0.
    seed_mismatches : int, optional
        Maximum mismatches allowed in the adapter seed used by ILLUMINACLIP.
        Default is 2.
    palindrome_clip_threshold : int, optional
        Threshold used by paired-end palindrome mode in ILLUMINACLIP.
        Default is 30.
    simple_clip_threshold : int, optional
        Threshold used by simple adapter clipping in ILLUMINACLIP.
        Default is 10.
    threads : int or {"auto"}, optional
        Number of Trimmomatic threads to request. Use 0 or "auto" to let
        Trimmomatic choose automatically. Default is 1.

    Returns
    -------
    CasavaOneEightSingleLanePerSampleDirFmt
        Trimmed paired-end sequences (both reads survived trimming).
    CasavaOneEightSingleLanePerSampleDirFmt
        Trimmed unpaired forward reads (reverse read failed trimming).
    CasavaOneEightSingleLanePerSampleDirFmt
        Trimmed unpaired reverse reads (forward read failed trimming).

    """
    adapter_path = files("q2_trimmomatic.bin.adapters").joinpath(adapter_file)
    executable_path = files("q2_trimmomatic.bin").joinpath("trimmomatic-0.39.jar")
    thread_args = _build_thread_args(threads)

    tasks = [
        (fwd, rev)
        for _, fwd, rev in paired_sequences.manifest.view(pd.DataFrame).itertuples()
    ]
    partitions = _run_parallel(
        tasks,
        n_jobs=n_jobs,
        func=lambda task: _trim_paired_partition(
            executable_path=executable_path,
            adapter_path=adapter_path,
            thread_args=thread_args,
            fwd=task[0],
            rev=task[1],
            head_crop=head_crop,
            seed_mismatches=seed_mismatches,
            palindrome_clip_threshold=palindrome_clip_threshold,
            simple_clip_threshold=simple_clip_threshold,
            leading=leading,
            trailing=trailing,
            sliding_window_size=sliding_window_size,
            sliding_window_quality=sliding_window_quality,
            crop=crop,
            min_length=min_length,
        ),
    )

    paired_end_trimmed = CasavaOneEightSingleLanePerSampleDirFmt()
    unpaired_fwd = CasavaOneEightSingleLanePerSampleDirFmt()
    unpaired_rev = CasavaOneEightSingleLanePerSampleDirFmt()
    for paired_partition, fwd_partition, rev_partition in partitions:
        _collate_partition(paired_partition, paired_end_trimmed)
        _collate_partition(fwd_partition, unpaired_fwd)
        _collate_partition(rev_partition, unpaired_rev)

    return paired_end_trimmed, unpaired_fwd, unpaired_rev


def trim_single(
    sequences: SingleLanePerSampleSingleEndFastqDirFmt,
    adapter_file: str = "TruSeq3-SE.fa",
    n_jobs: int = 1,
    leading: int = 3,
    trailing: int = 3,
    sliding_window_size: int = 4,
    sliding_window_quality: int = 15,
    min_length: int = 100,
    head_crop: int = 0,
    crop: int = 0,
    seed_mismatches: int = 2,
    simple_clip_threshold: int = 10,
    threads: int | str = 1,
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """Trim single-end reads using Trimmomatic.

    Parameters
    ----------
    sequences : SingleLanePerSampleSingleEndFastqDirFmt
        Single-end FASTQ sequence data to trim.
    adapter_file : str, optional
        Name of the bundled adapter FASTA file to use for ILLUMINACLIP.
        Default is "TruSeq3-SE.fa".
    n_jobs : int, optional
        Number of samples to process concurrently. Each concurrent job
        launches a separate Trimmomatic process. Default is 1.
    leading : int, optional
        Minimum quality required to keep a base at the start of a read (LEADING).
        Default is 3.
    trailing : int, optional
        Minimum quality required to keep a base at the end of a read (TRAILING).
        Default is 3.
    sliding_window_size : int, optional
        Number of bases in the sliding window for quality averaging (SLIDINGWINDOW).
        Default is 4.
    sliding_window_quality : int, optional
        Minimum average quality required in the sliding window (SLIDINGWINDOW).
        Default is 15.
    min_length : int, optional
        Reads shorter than this length after trimming will be discarded (MINLEN).
        Default is 100.
    head_crop : int, optional
        Number of bases to remove from the start of each read (HEADCROP).
        Applied before adapter clipping. Set to 0 to disable. Default is 0.
    crop : int, optional
        Cut reads to this length by removing bases from the end (CROP).
        Applied after quality trimming but before MINLEN. Set to 0 to disable.
        Default is 0.
    seed_mismatches : int, optional
        Maximum mismatches allowed in the adapter seed used by ILLUMINACLIP.
        Default is 2.
    simple_clip_threshold : int, optional
        Threshold used by simple adapter clipping in ILLUMINACLIP.
        Default is 10.
    threads : int or {"auto"}, optional
        Number of Trimmomatic threads to request. Use 0 or "auto" to let
        Trimmomatic choose automatically. Default is 1.

    Returns
    -------
    CasavaOneEightSingleLanePerSampleDirFmt
        Trimmed single-end sequences.

    """
    adapter_path = files("q2_trimmomatic.bin.adapters").joinpath(adapter_file)
    executable_path = files("q2_trimmomatic.bin").joinpath("trimmomatic-0.39.jar")
    thread_args = _build_thread_args(threads)

    tasks = [filepath for _, filepath in sequences.manifest.view(pd.DataFrame).itertuples()]
    partitions = _run_parallel(
        tasks,
        n_jobs=n_jobs,
        func=lambda filepath: _trim_single_partition(
            executable_path=executable_path,
            adapter_path=adapter_path,
            thread_args=thread_args,
            filepath=filepath,
            head_crop=head_crop,
            seed_mismatches=seed_mismatches,
            simple_clip_threshold=simple_clip_threshold,
            leading=leading,
            trailing=trailing,
            sliding_window_size=sliding_window_size,
            sliding_window_quality=sliding_window_quality,
            crop=crop,
            min_length=min_length,
        ),
    )

    trimmed = CasavaOneEightSingleLanePerSampleDirFmt()
    for partition in partitions:
        _collate_partition(partition, trimmed)

    return trimmed


def _trim_paired_partition(
    executable_path,
    adapter_path,
    thread_args,
    fwd,
    rev,
    head_crop,
    seed_mismatches,
    palindrome_clip_threshold,
    simple_clip_threshold,
    leading,
    trailing,
    sliding_window_size,
    sliding_window_quality,
    crop,
    min_length,
):
    """Run Trimmomatic on a single paired-end sample pair."""
    paired_end_trimmed = CasavaOneEightSingleLanePerSampleDirFmt()
    unpaired_fwd = CasavaOneEightSingleLanePerSampleDirFmt()
    unpaired_rev = CasavaOneEightSingleLanePerSampleDirFmt()

    cmd = [
        "java",
        "-jar",
        str(executable_path),
        "PE",
        *thread_args,
        fwd,
        rev,
        str(paired_end_trimmed.path / os.path.basename(fwd)),
        str(unpaired_fwd.path / os.path.basename(fwd)),
        str(paired_end_trimmed.path / os.path.basename(rev)),
        str(unpaired_rev.path / os.path.basename(rev)),
    ]
    if head_crop > 0:
        cmd.append(f"HEADCROP:{head_crop}")
    cmd += [
        f"ILLUMINACLIP:{adapter_path}:{seed_mismatches}:{palindrome_clip_threshold}:{simple_clip_threshold}",
        f"LEADING:{leading}",
        f"TRAILING:{trailing}",
        f"SLIDINGWINDOW:{sliding_window_size}:{sliding_window_quality}",
    ]
    if crop > 0:
        cmd.append(f"CROP:{crop}")
    cmd.append(f"MINLEN:{min_length}")

    subprocess.run(cmd, check=True)
    return paired_end_trimmed, unpaired_fwd, unpaired_rev


def _trim_single_partition(
    executable_path,
    adapter_path,
    thread_args,
    filepath,
    head_crop,
    seed_mismatches,
    simple_clip_threshold,
    leading,
    trailing,
    sliding_window_size,
    sliding_window_quality,
    crop,
    min_length,
):
    """Run Trimmomatic on a single single-end sample."""
    trimmed = CasavaOneEightSingleLanePerSampleDirFmt()

    cmd = [
        "java",
        "-jar",
        str(executable_path),
        "SE",
        *thread_args,
        filepath,
        str(trimmed.path / os.path.basename(filepath)),
    ]
    if head_crop > 0:
        cmd.append(f"HEADCROP:{head_crop}")
    cmd += [
        f"ILLUMINACLIP:{adapter_path}:{seed_mismatches}:{simple_clip_threshold}",
        f"LEADING:{leading}",
        f"TRAILING:{trailing}",
        f"SLIDINGWINDOW:{sliding_window_size}:{sliding_window_quality}",
    ]
    if crop > 0:
        cmd.append(f"CROP:{crop}")
    cmd.append(f"MINLEN:{min_length}")

    subprocess.run(cmd, check=True)
    return trimmed


def _run_parallel(tasks, n_jobs, func):
    """Run work serially or with a thread pool while preserving task order."""
    if n_jobs == 1:
        return [func(task) for task in tasks]

    max_workers = min(n_jobs, len(tasks)) if tasks else 1
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        return list(executor.map(func, tasks))


def _collate_partition(source, destination):
    """Copy all files from a partition output directory into the destination."""
    for fp in source.path.iterdir():
        shutil.copy2(fp, destination.path / fp.name)


def _build_thread_args(threads):
    """Translate QIIME ``Threads`` values to Trimmomatic CLI arguments."""
    if threads in (None, 0, "auto"):
        return []

    if isinstance(threads, str) and threads.isdigit():
        threads = int(threads)

    return ["-threads", str(threads)]
