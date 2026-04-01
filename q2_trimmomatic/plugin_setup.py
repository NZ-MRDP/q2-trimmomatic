"""QIIME 2 plugin setup for Trimmomatic."""

from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from qiime2.plugin import Collection, Choices, Int, List, Plugin, Range, Str, Threads

import q2_trimmomatic

from . import __version__

_ADAPTER_FILES_PE = [
    "NexteraPE-PE.fa",
    "TruSeq2-PE.fa",
    "TruSeq3-PE-2.fa",
    "TruSeq3-PE.fa",
]

_ADAPTER_FILES_SE = [
    "TruSeq2-SE.fa",
    "TruSeq3-SE.fa",
]

_ADAPTER_FILE_DESCRIPTION = (
    "Bundled adapter FASTA file to use for adapter clipping (ILLUMINACLIP). "
    "Choose the file matching your library preparation kit."
)

_SHARED_PARAMETER_DESCRIPTIONS = {
    "threads": (
        "Number of CPU threads to use per Trimmomatic process. Use 0 or auto "
        "to let Trimmomatic choose based on available resources."
    ),
    "leading": (
        "Minimum quality required to keep a base at the start of a read (LEADING). "
        "Bases below this threshold are removed."
    ),
    "trailing": (
        "Minimum quality required to keep a base at the end of a read (TRAILING). "
        "Bases below this threshold are removed."
    ),
    "sliding_window_size": ("Number of bases in the sliding window for average quality calculation (SLIDINGWINDOW)."),
    "sliding_window_quality": (
        "Minimum average quality required within the sliding window (SLIDINGWINDOW). "
        "The read is clipped once the average falls below this value."
    ),
    "min_length": "Reads shorter than this length after trimming will be discarded (MINLEN).",
    "head_crop": (
        "Number of bases to remove from the start of each read (HEADCROP). "
        "Useful for removing primer sequences. Set to 0 to disable."
    ),
    "crop": "Cut reads to this length by removing bases from the end (CROP). Set to 0 to disable.",
    "seed_mismatches": (
        "Maximum mismatches allowed in the adapter seed during ILLUMINACLIP. Higher values are more permissive."
    ),
    "simple_clip_threshold": "Threshold for simple adapter clipping in ILLUMINACLIP.",
}

_SHARED_PARAMETERS = {
    "threads": Threads,
    "leading": Int % Range(0, None),
    "trailing": Int % Range(0, None),
    "sliding_window_size": Int % Range(1, None),
    "sliding_window_quality": Int % Range(0, None),
    "min_length": Int % Range(1, None),
    "head_crop": Int % Range(0, None),
    "crop": Int % Range(0, None),
    "seed_mismatches": Int % Range(0, None),
    "simple_clip_threshold": Int % Range(0, None),
}

_PAIRED_ONLY_PARAMETER_DESCRIPTIONS = {
    "palindrome_clip_threshold": (
        "Threshold for palindrome adapter clipping in paired-end ILLUMINACLIP mode."
    ),
}

_PAIRED_ONLY_PARAMETERS = {
    "palindrome_clip_threshold": Int % Range(0, None),
}

plugin = Plugin(
    name="trimmomatic",
    version=__version__,
    description="QIIME 2 plugin for adapter trimming and quality control with Trimmomatic.",
    website="https://github.com/NZ-MRDP/q2-trimmomatic",
    package="q2_trimmomatic",
    user_support_text="Please post questions to the QIIME 2 forum at https://forum.qiime2.org.",
)

plugin.methods.register_function(
    function=q2_trimmomatic.trim_paired,
    name="Remove adapters from paired sequences and quality trim",
    description="Remove adapters from paired-end sequences and perform quality trimming using Trimmomatic.",
    inputs={
        "paired_sequences": SampleData[PairedEndSequencesWithQuality],
    },
    parameters={
        "adapter_file": Str % Choices(_ADAPTER_FILES_PE),
        **_SHARED_PARAMETERS,
        **_PAIRED_ONLY_PARAMETERS,
    },
    outputs=[
        ("paired_end_trimmed", SampleData[PairedEndSequencesWithQuality]),
        ("unpaired_fwd", SampleData[SequencesWithQuality]),
        ("unpaired_rev", SampleData[SequencesWithQuality]),
    ],
    input_descriptions={
        "paired_sequences": "Illumina paired-end sequence data.",
    },
    parameter_descriptions={
        "adapter_file": _ADAPTER_FILE_DESCRIPTION,
        **_SHARED_PARAMETER_DESCRIPTIONS,
        **_PAIRED_ONLY_PARAMETER_DESCRIPTIONS,
    },
    output_descriptions={
        "paired_end_trimmed": "Trimmed paired-end reads where both mates passed filtering.",
        "unpaired_fwd": "Trimmed forward reads whose reverse mate failed filtering.",
        "unpaired_rev": "Trimmed reverse reads whose forward mate failed filtering.",
    },
)

plugin.methods.register_function(
    function=q2_trimmomatic.trim_single,
    name="Remove adapters from single-end sequences and quality trim",
    description="Remove adapters from single-end sequences and perform quality trimming using Trimmomatic.",
    inputs={
        "sequences": SampleData[SequencesWithQuality],
    },
    parameters={
        "adapter_file": Str % Choices(_ADAPTER_FILES_SE),
        **_SHARED_PARAMETERS,
    },
    outputs=[
        ("trimmed", SampleData[SequencesWithQuality]),
    ],
    input_descriptions={
        "sequences": "Illumina single-end sequence data.",
    },
    parameter_descriptions={
        "adapter_file": _ADAPTER_FILE_DESCRIPTION,
        **_SHARED_PARAMETER_DESCRIPTIONS,
    },
    output_descriptions={
        "trimmed": "Trimmed single-end reads that passed the requested adapter and quality filters.",
    },
)

plugin.methods.register_function(
    function=q2_trimmomatic.split_paired_samples,
    name="Split paired-end samples into partitions",
    description=(
        "Split paired-end sequence data into a collection of smaller "
        "artifacts for use in the plugin's parallel trimming pipeline."
    ),
    inputs={
        "paired_sequences": SampleData[PairedEndSequencesWithQuality],
    },
    parameters={
        "num_partitions": Int % Range(1, None),
    },
    outputs=[
        ("partitions", Collection[SampleData[PairedEndSequencesWithQuality]]),
    ],
    input_descriptions={
        "paired_sequences": "Illumina paired-end sequence data.",
    },
    parameter_descriptions={
        "num_partitions": (
            "Number of partitions to split samples into. Defaults to one "
            "partition per sample."
        ),
    },
    output_descriptions={
        "partitions": "A collection of smaller paired-end sequence artifacts.",
    },
)

plugin.methods.register_function(
    function=q2_trimmomatic.split_single_samples,
    name="Split single-end samples into partitions",
    description=(
        "Split single-end sequence data into a collection of smaller "
        "artifacts for use in the plugin's parallel trimming pipeline."
    ),
    inputs={
        "sequences": SampleData[SequencesWithQuality],
    },
    parameters={
        "num_partitions": Int % Range(1, None),
    },
    outputs=[
        ("partitions", Collection[SampleData[SequencesWithQuality]]),
    ],
    input_descriptions={
        "sequences": "Illumina single-end sequence data.",
    },
    parameter_descriptions={
        "num_partitions": (
            "Number of partitions to split samples into. Defaults to one "
            "partition per sample."
        ),
    },
    output_descriptions={
        "partitions": "A collection of smaller single-end sequence artifacts.",
    },
)

plugin.methods.register_function(
    function=q2_trimmomatic.collate_trimmed_paired,
    name="Collate trimmed paired-end sequence partitions",
    description=(
        "Collate a collection of trimmed paired-end sequence artifacts "
        "(produced by a parallel pipeline) into a single artifact."
    ),
    inputs={
        "trimmed": List[SampleData[PairedEndSequencesWithQuality]],
    },
    parameters={},
    outputs=[
        ("collated", SampleData[PairedEndSequencesWithQuality]),
    ],
    input_descriptions={
        "trimmed": "A list of trimmed paired-end sequence artifacts to collate.",
    },
    parameter_descriptions={},
    output_descriptions={
        "collated": "The collated trimmed paired-end sequences.",
    },
)

plugin.methods.register_function(
    function=q2_trimmomatic.collate_trimmed_single,
    name="Collate trimmed single-end sequence partitions",
    description=(
        "Collate a collection of trimmed single-end sequence artifacts "
        "(produced by a parallel pipeline) into a single artifact."
    ),
    inputs={
        "trimmed": List[SampleData[SequencesWithQuality]],
    },
    parameters={},
    outputs=[
        ("collated", SampleData[SequencesWithQuality]),
    ],
    input_descriptions={
        "trimmed": "A list of trimmed single-end sequence artifacts to collate.",
    },
    parameter_descriptions={},
    output_descriptions={
        "collated": "The collated trimmed single-end sequences.",
    },
)

plugin.pipelines.register_function(
    function=q2_trimmomatic.trim_paired_parallel,
    name="Parallel adapter trimming of paired-end sequences",
    description=(
        "Partition samples, trim each partition in parallel using "
        "trimmomatic.trim_paired, then collate the results. When a parallel "
        "execution backend is configured (e.g. parsl), partitions are "
        "processed concurrently."
    ),
    inputs={
        "paired_sequences": SampleData[PairedEndSequencesWithQuality],
    },
    parameters={
        "num_partitions": Int % Range(1, None),
        "adapter_file": Str % Choices(_ADAPTER_FILES_PE),
        **_SHARED_PARAMETERS,
        **_PAIRED_ONLY_PARAMETERS,
    },
    outputs=[
        ("paired_end_trimmed", SampleData[PairedEndSequencesWithQuality]),
        ("unpaired_fwd", SampleData[SequencesWithQuality]),
        ("unpaired_rev", SampleData[SequencesWithQuality]),
    ],
    input_descriptions={
        "paired_sequences": "Illumina paired-end sequence data.",
    },
    parameter_descriptions={
        "num_partitions": (
            "Number of partitions to split samples into for parallel "
            "processing. Defaults to one partition per sample."
        ),
        "adapter_file": _ADAPTER_FILE_DESCRIPTION,
        **_SHARED_PARAMETER_DESCRIPTIONS,
        **_PAIRED_ONLY_PARAMETER_DESCRIPTIONS,
    },
    output_descriptions={
        "paired_end_trimmed": "Trimmed paired-end reads where both mates passed filtering.",
        "unpaired_fwd": "Trimmed forward reads whose reverse mate failed filtering.",
        "unpaired_rev": "Trimmed reverse reads whose forward mate failed filtering.",
    },
)

plugin.pipelines.register_function(
    function=q2_trimmomatic.trim_single_parallel,
    name="Parallel adapter trimming of single-end sequences",
    description=(
        "Partition samples, trim each partition in parallel using "
        "trimmomatic.trim_single, then collate the results. When a parallel "
        "execution backend is configured (e.g. parsl), partitions are "
        "processed concurrently."
    ),
    inputs={
        "sequences": SampleData[SequencesWithQuality],
    },
    parameters={
        "num_partitions": Int % Range(1, None),
        "adapter_file": Str % Choices(_ADAPTER_FILES_SE),
        **_SHARED_PARAMETERS,
    },
    outputs=[
        ("trimmed", SampleData[SequencesWithQuality]),
    ],
    input_descriptions={
        "sequences": "Illumina single-end sequence data.",
    },
    parameter_descriptions={
        "num_partitions": (
            "Number of partitions to split samples into for parallel "
            "processing. Defaults to one partition per sample."
        ),
        "adapter_file": _ADAPTER_FILE_DESCRIPTION,
        **_SHARED_PARAMETER_DESCRIPTIONS,
    },
    output_descriptions={
        "trimmed": "Trimmed single-end reads that passed the requested adapter and quality filters.",
    },
)
