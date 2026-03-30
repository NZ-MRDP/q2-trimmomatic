"""Trimmomatic adapter trimming and quality control functions."""

import os
import subprocess
from importlib.resources import files
from typing import Tuple

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)


def trim_paired(
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
) -> Tuple[
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
]:
    """Trim paired-end reads using Trimmomatic.

    Parameters
    ----------
    paired_sequences : SingleLanePerSamplePairedEndFastqDirFmt
        Paired-end FASTQ sequence data to trim.
    adapter_file : str, optional
        Name of the bundled adapter FASTA file to use for ILLUMINACLIP.
        Default is "NexteraPE-PE.fa".
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

    Returns
    -------
    SingleLanePerSamplePairedEndFastqDirFmt
        Trimmed paired-end sequences (both reads survived trimming).
    SingleLanePerSampleSingleEndFastqDirFmt
        Trimmed unpaired forward reads (reverse read failed trimming).
    SingleLanePerSampleSingleEndFastqDirFmt
        Trimmed unpaired reverse reads (forward read failed trimming).

    """
    adapter_path = files("q2_trimmomatic.bin.adapters").joinpath(adapter_file)
    executable_path = files("q2_trimmomatic.bin").joinpath("trimmomatic-0.39.jar")

    paired_end_trimmed = SingleLanePerSamplePairedEndFastqDirFmt()
    unpaired_fwd = SingleLanePerSampleSingleEndFastqDirFmt()
    unpaired_rev = SingleLanePerSampleSingleEndFastqDirFmt()

    df = paired_sequences.manifest.view(pd.DataFrame)
    for _, fwd, rev in df.itertuples():
        trimmed_paired_fwd = os.path.join(str(paired_end_trimmed), os.path.basename(fwd))
        trimmed_paired_rev = os.path.join(str(paired_end_trimmed), os.path.basename(rev))
        trimmed_unpaired_fwd = os.path.join(str(unpaired_fwd), os.path.basename(fwd))
        trimmed_unpaired_rev = os.path.join(str(unpaired_rev), os.path.basename(rev))

        # Trimmomatic processes steps in order. HEADCROP runs first (before adapter
        # clipping), CROP runs after quality trimming but before MINLEN.
        cmd = [
            "java",
            "-jar",
            str(executable_path),
            "PE",
            fwd,
            rev,
            trimmed_paired_fwd,
            trimmed_unpaired_fwd,
            trimmed_paired_rev,
            trimmed_unpaired_rev,
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


def trim_single(
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
) -> SingleLanePerSampleSingleEndFastqDirFmt:
    """Trim single-end reads using Trimmomatic.

    Parameters
    ----------
    sequences : SingleLanePerSampleSingleEndFastqDirFmt
        Single-end FASTQ sequence data to trim.
    adapter_file : str, optional
        Name of the bundled adapter FASTA file to use for ILLUMINACLIP.
        Default is "TruSeq3-SE.fa".
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

    Returns
    -------
    SingleLanePerSampleSingleEndFastqDirFmt
        Trimmed single-end sequences.

    """
    adapter_path = files("q2_trimmomatic.bin.adapters").joinpath(adapter_file)
    executable_path = files("q2_trimmomatic.bin").joinpath("trimmomatic-0.39.jar")

    trimmed = SingleLanePerSampleSingleEndFastqDirFmt()

    df = sequences.manifest.view(pd.DataFrame)
    for _, filepath in df.itertuples():
        trimmed_path = os.path.join(str(trimmed), os.path.basename(filepath))

        # Trimmomatic SE ILLUMINACLIP takes seedMismatches:simpleClipThreshold
        # (no palindromeClipThreshold — that is PE-only).
        cmd = [
            "java",
            "-jar",
            str(executable_path),
            "SE",
            filepath,
            trimmed_path,
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


trim_paired.__annotations__["return"] = (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
