import os
import subprocess
from importlib import resources

import pandas as pd
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)

from q2_trimmomatic import bin
from q2_trimmomatic.bin import adapters


def trim_paired(paired_sequences: SingleLanePerSamplePairedEndFastqDirFmt, min_length: int = 100) -> (
    CasavaOneEightSingleLanePerSampleDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt,
):  # type: ignore
    """
    Purpose:  trim paired end reads using trimmomatic

    Inputs:
    ----------
    pe1_fastq_gz:
        raw fastq.gz read files (read 1)

    pe2_fastq_gz:
        raw fastq.gz read files (read 2)

    sample_base_name:
        user-specified descriptor for the sample, which will be
         used to label all relevant output files and the file directory

     trimmomatic_path:
         filepath to trimmomatic

     adapter_path:
         filepath to adaptors

    Returns
    ----------
    /Trimmed reads/:
        directory containing output trimmed reads files
    trimmed_pe1:
        trimmed paired end 1 reads
    trimmed_pe2:
        trimmed paired end 2 reads
    trimmed_u1:
        trimmed unpaired (orphan) end 1 reads
    trimmed_u2:
        trimmed unpaired (orphan) end 2 reads

    External dependencies
    ----------
    Trimmomatic:
        http://www.usadellab.org/cms/?page=trimmomatic


    Notes
    ----------
    1. Unless this script and all of the input arguments are by default in your
         $PATH, you will need to include the PATH to each of these in your command
    2. This script has been written to run within a conda environment containing
         the necessary dependencies.

    Dependencies
        import sys
        import subprocess
        make_new_dir

    """
    #### Step 1:  Trim reads
    # new variables defined within this code block:
    # trimmed_pe1 = name of trimmed read 1 (paired, *.fastq.gz)
    # trimmed_pe2 = name of trimmed read 2 (paired, *.fastq.gz)
    # trimmed_u1 = name of trimmed read 1 (unpaired, *.fastq.gz)
    # trimmed_u2 = name of trimmed read 2 (unpaired, *.fastq.gz)
    # Trimmed_reads = folder of trimmed reads

    # temp_dir = "./Trimmed_reads/"
    # make_new_dir(temp_dir)

    # trimmed_pe1 = f'{temp_dir}{pe1_fastq_gz.split("/")[-1].replace(".fastq.gz", "_trimmed_pe1.fastq.gz")}'
    # trimmed_pe2 = f'{temp_dir}{pe2_fastq_gz.split("/")[-1].replace(".fastq.gz", "_trimmed_pe2.fastq.gz")}'
    # trimmed_u1 = f'{temp_dir}{pe1_fastq_gz.split("/")[-1].replace(".fastq.gz", "_trimmed_u1.fastq.gz")}'
    # trimmed_u2 = f'{temp_dir}{pe2_fastq_gz.split("/")[-1].replace(".fastq.gz", "_trimmed_u2.fastq.gz")}'

    # syscommand = f"java -jar {trimmomatic_path} PE {pe1_fastq_gz} {pe2_fastq_gz} {trimmed_pe1} {trimmed_u1} {trimmed_pe2} {trimmed_u2} ILLUMINACLIP:{adapter_path}NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100"
    # subprocess.call(syscommand, shell=True)
    paired_end_trimmed = CasavaOneEightSingleLanePerSampleDirFmt()
    unpaired_fwd = CasavaOneEightSingleLanePerSampleDirFmt()
    unpaired_rev = CasavaOneEightSingleLanePerSampleDirFmt()
    # TODO: make adapted file path into an input
    with resources.path(adapters, "NexteraPE-PE.fa") as adapter_path, resources.path(
        bin, "trimmomatic-0.39.jar"
    ) as executable_path:
        df = paired_sequences.manifest.view(pd.DataFrame)
        for _, fwd, rev in df.itertuples():
            trimmed_paired_fwd = os.path.join(str(paired_end_trimmed), os.path.basename(fwd))
            trimmed_paired_rev = os.path.join(str(paired_end_trimmed), os.path.basename(rev))
            trimmed_unpaired_fwd = os.path.join(str(unpaired_fwd), os.path.basename(fwd))
            trimmed_unpaired_rev = os.path.join(str(unpaired_rev), os.path.basename(rev))
            subprocess.run(
                [
                    "java",
                    "-jar",
                    executable_path,
                    "PE",
                    fwd,
                    rev,
                    trimmed_paired_fwd,
                    trimmed_unpaired_fwd,
                    trimmed_paired_rev,
                    trimmed_unpaired_rev,
                    f"ILLUMINACLIP:{adapter_path}:2:30:10",
                    "LEADING:3",
                    "TRAILING:3",
                    "SLIDINGWINDOW:4:15",
                    f"MINLEN:{min_length}",
                ]
            )
    return paired_end_trimmed, unpaired_fwd, unpaired_rev
