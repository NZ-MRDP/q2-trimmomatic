"""QIIME 2 plugin setup for Trimmomatic."""

from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from qiime2.plugin import Choices, Int, Plugin, Range, Str

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
}

_SHARED_PARAMETERS = {
    "leading": Int % Range(0, None),
    "trailing": Int % Range(0, None),
    "sliding_window_size": Int % Range(1, None),
    "sliding_window_quality": Int % Range(0, None),
    "min_length": Int % Range(1, None),
    "head_crop": Int % Range(0, None),
    "crop": Int % Range(0, None),
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
    },
    output_descriptions={
        "paired_end_trimmed": "Trimmed paired-end sequence data.",
        "unpaired_fwd": "Trimmed unpaired forward reads.",
        "unpaired_rev": "Trimmed unpaired reverse reads.",
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
        "trimmed": "Trimmed single-end sequence data.",
    },
)
