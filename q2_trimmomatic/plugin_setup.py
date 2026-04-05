"""QIIME 2 plugin setup for Trimmomatic."""

from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from qiime2.plugin import Collection, Choices, Int, Plugin, Range, Str, Threads

from . import __version__
from ._pipelines import (
    _combine_paired_sequences,
    _combine_single_sequences,
    _split_paired_samples,
    _split_single_samples,
    trim_paired,
    trim_single,
)
from ._trimmomatic import _trim_paired_worker, _trim_single_worker

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

_PIPELINE_PARAMETER_DESCRIPTIONS = {
    "num_partitions": (
        "Number of internal partitions to create for split/apply/combine "
        "execution. Defaults to one partition per sample."
    ),
    "threads": (
        "Number of CPU threads to use per Trimmomatic worker process. Use 0 "
        "or auto to let Trimmomatic choose based on available resources."
    ),
    "leading": (
        "Minimum quality required to keep a base at the start of a read "
        "(LEADING). Bases below this threshold are removed."
    ),
    "trailing": (
        "Minimum quality required to keep a base at the end of a read "
        "(TRAILING). Bases below this threshold are removed."
    ),
    "sliding_window_size": (
        "Number of bases in the sliding window for average quality "
        "calculation (SLIDINGWINDOW)."
    ),
    "sliding_window_quality": (
        "Minimum average quality required within the sliding window "
        "(SLIDINGWINDOW). The read is clipped once the average falls below "
        "this value."
    ),
    "min_length": (
        "Reads shorter than this length after trimming will be discarded "
        "(MINLEN)."
    ),
    "head_crop": (
        "Number of bases to remove from the start of each read (HEADCROP). "
        "Useful for removing primer sequences. Set to 0 to disable."
    ),
    "crop": (
        "Cut reads to this length by removing bases from the end (CROP). "
        "Set to 0 to disable."
    ),
    "seed_mismatches": (
        "Maximum mismatches allowed in the adapter seed during "
        "ILLUMINACLIP. Higher values are more permissive."
    ),
    "simple_clip_threshold": (
        "Threshold for simple adapter clipping in ILLUMINACLIP."
    ),
}

_PIPELINE_PARAMETERS = {
    "num_partitions": Int % Range(1, None),
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
        "Threshold for palindrome adapter clipping in paired-end "
        "ILLUMINACLIP mode."
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
    function=_split_paired_samples,
    name="Internal: split paired-end samples",
    description=(
        "Internal helper for the paired-end trimming pipeline. Splits one "
        "paired-end artifact into a collection of partition artifacts."
    ),
    inputs={
        "paired_sequences": SampleData[PairedEndSequencesWithQuality],
    },
    parameters={
        "num_partitions": Int % Range(1, None),
    },
    outputs=[
        ("pieces", Collection[SampleData[PairedEndSequencesWithQuality]]),
    ],
    input_descriptions={
        "paired_sequences": "Illumina paired-end sequence data.",
    },
    parameter_descriptions={
        "num_partitions": _PIPELINE_PARAMETER_DESCRIPTIONS["num_partitions"],
    },
    output_descriptions={
        "pieces": "Internal collection of paired-end partitions.",
    },
)

plugin.methods.register_function(
    function=_split_single_samples,
    name="Internal: split single-end samples",
    description=(
        "Internal helper for the single-end trimming pipeline. Splits one "
        "single-end artifact into a collection of partition artifacts."
    ),
    inputs={
        "sequences": SampleData[SequencesWithQuality],
    },
    parameters={
        "num_partitions": Int % Range(1, None),
    },
    outputs=[
        ("pieces", Collection[SampleData[SequencesWithQuality]]),
    ],
    input_descriptions={
        "sequences": "Illumina single-end sequence data.",
    },
    parameter_descriptions={
        "num_partitions": _PIPELINE_PARAMETER_DESCRIPTIONS["num_partitions"],
    },
    output_descriptions={
        "pieces": "Internal collection of single-end partitions.",
    },
)

plugin.methods.register_function(
    function=_trim_paired_worker,
    name="Internal: trim paired-end partition",
    description=(
        "Internal worker action for the paired-end trimming pipeline. Trims "
        "one paired-end partition artifact."
    ),
    inputs={
        "paired_sequences": SampleData[PairedEndSequencesWithQuality],
    },
    parameters={
        "adapter_file": Str % Choices(_ADAPTER_FILES_PE),
        "threads": Threads,
        "leading": Int % Range(0, None),
        "trailing": Int % Range(0, None),
        "sliding_window_size": Int % Range(1, None),
        "sliding_window_quality": Int % Range(0, None),
        "min_length": Int % Range(1, None),
        "head_crop": Int % Range(0, None),
        "crop": Int % Range(0, None),
        "seed_mismatches": Int % Range(0, None),
        "palindrome_clip_threshold": Int % Range(0, None),
        "simple_clip_threshold": Int % Range(0, None),
    },
    outputs=[
        ("paired_end_trimmed", SampleData[PairedEndSequencesWithQuality]),
        ("unpaired_fwd", SampleData[SequencesWithQuality]),
        ("unpaired_rev", SampleData[SequencesWithQuality]),
    ],
    input_descriptions={
        "paired_sequences": "One internal paired-end partition artifact.",
    },
    parameter_descriptions={
        "adapter_file": _ADAPTER_FILE_DESCRIPTION,
        "threads": _PIPELINE_PARAMETER_DESCRIPTIONS["threads"],
        "leading": _PIPELINE_PARAMETER_DESCRIPTIONS["leading"],
        "trailing": _PIPELINE_PARAMETER_DESCRIPTIONS["trailing"],
        "sliding_window_size": _PIPELINE_PARAMETER_DESCRIPTIONS["sliding_window_size"],
        "sliding_window_quality": _PIPELINE_PARAMETER_DESCRIPTIONS["sliding_window_quality"],
        "min_length": _PIPELINE_PARAMETER_DESCRIPTIONS["min_length"],
        "head_crop": _PIPELINE_PARAMETER_DESCRIPTIONS["head_crop"],
        "crop": _PIPELINE_PARAMETER_DESCRIPTIONS["crop"],
        "seed_mismatches": _PIPELINE_PARAMETER_DESCRIPTIONS["seed_mismatches"],
        "palindrome_clip_threshold": _PAIRED_ONLY_PARAMETER_DESCRIPTIONS["palindrome_clip_threshold"],
        "simple_clip_threshold": _PIPELINE_PARAMETER_DESCRIPTIONS["simple_clip_threshold"],
    },
    output_descriptions={
        "paired_end_trimmed": "Trimmed paired-end reads where both mates passed filtering.",
        "unpaired_fwd": "Trimmed forward reads whose reverse mate failed filtering.",
        "unpaired_rev": "Trimmed reverse reads whose forward mate failed filtering.",
    },
)

plugin.methods.register_function(
    function=_trim_single_worker,
    name="Internal: trim single-end partition",
    description=(
        "Internal worker action for the single-end trimming pipeline. Trims "
        "one single-end partition artifact."
    ),
    inputs={
        "sequences": SampleData[SequencesWithQuality],
    },
    parameters={
        "adapter_file": Str % Choices(_ADAPTER_FILES_SE),
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
    },
    outputs=[
        ("trimmed", SampleData[SequencesWithQuality]),
    ],
    input_descriptions={
        "sequences": "One internal single-end partition artifact.",
    },
    parameter_descriptions={
        "adapter_file": _ADAPTER_FILE_DESCRIPTION,
        "threads": _PIPELINE_PARAMETER_DESCRIPTIONS["threads"],
        "leading": _PIPELINE_PARAMETER_DESCRIPTIONS["leading"],
        "trailing": _PIPELINE_PARAMETER_DESCRIPTIONS["trailing"],
        "sliding_window_size": _PIPELINE_PARAMETER_DESCRIPTIONS["sliding_window_size"],
        "sliding_window_quality": _PIPELINE_PARAMETER_DESCRIPTIONS["sliding_window_quality"],
        "min_length": _PIPELINE_PARAMETER_DESCRIPTIONS["min_length"],
        "head_crop": _PIPELINE_PARAMETER_DESCRIPTIONS["head_crop"],
        "crop": _PIPELINE_PARAMETER_DESCRIPTIONS["crop"],
        "seed_mismatches": _PIPELINE_PARAMETER_DESCRIPTIONS["seed_mismatches"],
        "simple_clip_threshold": _PIPELINE_PARAMETER_DESCRIPTIONS["simple_clip_threshold"],
    },
    output_descriptions={
        "trimmed": "Trimmed single-end reads that passed the requested adapter and quality filters.",
    },
)

plugin.methods.register_function(
    function=_combine_paired_sequences,
    name="Internal: combine paired-end partitions",
    description=(
        "Internal helper for the paired-end trimming pipeline. Combines a "
        "collection of paired-end partition outputs into one artifact."
    ),
    inputs={
        "trimmed": Collection[SampleData[PairedEndSequencesWithQuality]],
    },
    parameters={},
    outputs=[
        ("combined", SampleData[PairedEndSequencesWithQuality]),
    ],
    input_descriptions={
        "trimmed": "Internal collection of trimmed paired-end partitions.",
    },
    parameter_descriptions={},
    output_descriptions={
        "combined": "Collated paired-end trimming results.",
    },
)

plugin.methods.register_function(
    function=_combine_single_sequences,
    name="Internal: combine single-end partitions",
    description=(
        "Internal helper for the trimming pipelines. Combines a collection "
        "of single-end partition outputs into one artifact."
    ),
    inputs={
        "trimmed": Collection[SampleData[SequencesWithQuality]],
    },
    parameters={},
    outputs=[
        ("combined", SampleData[SequencesWithQuality]),
    ],
    input_descriptions={
        "trimmed": "Internal collection of trimmed single-end partitions.",
    },
    parameter_descriptions={},
    output_descriptions={
        "combined": "Collated single-end trimming results.",
    },
)

plugin.pipelines.register_function(
    function=trim_paired,
    name="Remove adapters from paired sequences and quality trim",
    description=(
        "Remove adapters from paired-end sequences and perform quality "
        "trimming using Trimmomatic. This user-facing Pipeline uses an "
        "internal split/apply/combine implementation so the public interface "
        "remains normal QIIME 2 artifacts."
    ),
    inputs={
        "paired_sequences": SampleData[PairedEndSequencesWithQuality],
    },
    parameters={
        "num_partitions": Int % Range(1, None),
        "adapter_file": Str % Choices(_ADAPTER_FILES_PE),
        "threads": Threads,
        "leading": Int % Range(0, None),
        "trailing": Int % Range(0, None),
        "sliding_window_size": Int % Range(1, None),
        "sliding_window_quality": Int % Range(0, None),
        "min_length": Int % Range(1, None),
        "head_crop": Int % Range(0, None),
        "crop": Int % Range(0, None),
        "seed_mismatches": Int % Range(0, None),
        "palindrome_clip_threshold": Int % Range(0, None),
        "simple_clip_threshold": Int % Range(0, None),
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
        "num_partitions": _PIPELINE_PARAMETER_DESCRIPTIONS["num_partitions"],
        "adapter_file": _ADAPTER_FILE_DESCRIPTION,
        "threads": _PIPELINE_PARAMETER_DESCRIPTIONS["threads"],
        "leading": _PIPELINE_PARAMETER_DESCRIPTIONS["leading"],
        "trailing": _PIPELINE_PARAMETER_DESCRIPTIONS["trailing"],
        "sliding_window_size": _PIPELINE_PARAMETER_DESCRIPTIONS["sliding_window_size"],
        "sliding_window_quality": _PIPELINE_PARAMETER_DESCRIPTIONS["sliding_window_quality"],
        "min_length": _PIPELINE_PARAMETER_DESCRIPTIONS["min_length"],
        "head_crop": _PIPELINE_PARAMETER_DESCRIPTIONS["head_crop"],
        "crop": _PIPELINE_PARAMETER_DESCRIPTIONS["crop"],
        "seed_mismatches": _PIPELINE_PARAMETER_DESCRIPTIONS["seed_mismatches"],
        "palindrome_clip_threshold": _PAIRED_ONLY_PARAMETER_DESCRIPTIONS["palindrome_clip_threshold"],
        "simple_clip_threshold": _PIPELINE_PARAMETER_DESCRIPTIONS["simple_clip_threshold"],
    },
    output_descriptions={
        "paired_end_trimmed": "Trimmed paired-end reads where both mates passed filtering.",
        "unpaired_fwd": "Trimmed forward reads whose reverse mate failed filtering.",
        "unpaired_rev": "Trimmed reverse reads whose forward mate failed filtering.",
    },
)

plugin.pipelines.register_function(
    function=trim_single,
    name="Remove adapters from single-end sequences and quality trim",
    description=(
        "Remove adapters from single-end sequences and perform quality "
        "trimming using Trimmomatic. This user-facing Pipeline uses an "
        "internal split/apply/combine implementation so the public interface "
        "remains normal QIIME 2 artifacts."
    ),
    inputs={
        "sequences": SampleData[SequencesWithQuality],
    },
    parameters={
        "num_partitions": Int % Range(1, None),
        "adapter_file": Str % Choices(_ADAPTER_FILES_SE),
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
    },
    outputs=[
        ("trimmed", SampleData[SequencesWithQuality]),
    ],
    input_descriptions={
        "sequences": "Illumina single-end sequence data.",
    },
    parameter_descriptions={
        "num_partitions": _PIPELINE_PARAMETER_DESCRIPTIONS["num_partitions"],
        "adapter_file": _ADAPTER_FILE_DESCRIPTION,
        "threads": _PIPELINE_PARAMETER_DESCRIPTIONS["threads"],
        "leading": _PIPELINE_PARAMETER_DESCRIPTIONS["leading"],
        "trailing": _PIPELINE_PARAMETER_DESCRIPTIONS["trailing"],
        "sliding_window_size": _PIPELINE_PARAMETER_DESCRIPTIONS["sliding_window_size"],
        "sliding_window_quality": _PIPELINE_PARAMETER_DESCRIPTIONS["sliding_window_quality"],
        "min_length": _PIPELINE_PARAMETER_DESCRIPTIONS["min_length"],
        "head_crop": _PIPELINE_PARAMETER_DESCRIPTIONS["head_crop"],
        "crop": _PIPELINE_PARAMETER_DESCRIPTIONS["crop"],
        "seed_mismatches": _PIPELINE_PARAMETER_DESCRIPTIONS["seed_mismatches"],
        "simple_clip_threshold": _PIPELINE_PARAMETER_DESCRIPTIONS["simple_clip_threshold"],
    },
    output_descriptions={
        "trimmed": "Trimmed single-end reads that passed the requested adapter and quality filters.",
    },
)
