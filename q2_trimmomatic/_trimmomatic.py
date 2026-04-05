"""Trimmomatic worker methods used by the public QIIME 2 pipelines."""

import os
import subprocess
from importlib.resources import files

import pandas as pd
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)


def _trim_paired_worker(
    paired_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
    adapter_file: str = "NexteraPE-PE.fa",
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
) -> (
    CasavaOneEightSingleLanePerSampleDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt,
):
    """Trim one paired-end partition using Trimmomatic."""
    adapter_path = files("q2_trimmomatic.bin.adapters").joinpath(adapter_file)
    executable_path = files("q2_trimmomatic.bin").joinpath("trimmomatic-0.39.jar")
    thread_args = _build_thread_args(threads)

    paired_end_trimmed = CasavaOneEightSingleLanePerSampleDirFmt()
    unpaired_fwd = CasavaOneEightSingleLanePerSampleDirFmt()
    unpaired_rev = CasavaOneEightSingleLanePerSampleDirFmt()

    df = paired_sequences.manifest.view(pd.DataFrame)
    for _, fwd, rev in df.itertuples():
        fwd = str(_resolve_manifest_path(paired_sequences.path, fwd))
        rev = str(_resolve_manifest_path(paired_sequences.path, rev))
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
            (
                "ILLUMINACLIP:"
                f"{adapter_path}:{seed_mismatches}:"
                f"{palindrome_clip_threshold}:{simple_clip_threshold}"
            ),
            f"LEADING:{leading}",
            f"TRAILING:{trailing}",
            f"SLIDINGWINDOW:{sliding_window_size}:{sliding_window_quality}",
        ]
        if crop > 0:
            cmd.append(f"CROP:{crop}")
        cmd.append(f"MINLEN:{min_length}")

        subprocess.run(cmd, check=True)

    return paired_end_trimmed, unpaired_fwd, unpaired_rev


def _trim_single_worker(
    sequences: SingleLanePerSampleSingleEndFastqDirFmt,
    adapter_file: str = "TruSeq3-SE.fa",
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
    """Trim one single-end partition using Trimmomatic."""
    adapter_path = files("q2_trimmomatic.bin.adapters").joinpath(adapter_file)
    executable_path = files("q2_trimmomatic.bin").joinpath("trimmomatic-0.39.jar")
    thread_args = _build_thread_args(threads)

    trimmed = CasavaOneEightSingleLanePerSampleDirFmt()

    df = sequences.manifest.view(pd.DataFrame)
    for _, filepath in df.itertuples():
        filepath = str(_resolve_manifest_path(sequences.path, filepath))
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


def _build_thread_args(threads):
    """Translate QIIME ``Threads`` values to Trimmomatic CLI arguments."""
    if threads in (None, 0, "auto"):
        return []

    if isinstance(threads, str) and threads.isdigit():
        threads = int(threads)

    return ["-threads", str(threads)]


def _resolve_manifest_path(base_dir, filepath):
    """Resolve relative manifest entries against the source artifact path."""
    path = os.fspath(filepath)
    if os.path.isabs(path):
        return path

    return os.path.join(os.fspath(base_dir), path)
