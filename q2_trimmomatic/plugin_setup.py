from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from qiime2.plugin import Int, Plugin, Range

import q2_trimmomatic

from . import __version__

plugin = Plugin(
    name="trimmomatic",
    version=__version__,
    description="QIIME2 plugin for trimmomatic",
    website="https://github.com/NZ-MRDP/q2-trimmomatic",
    package="q2_trimmomatic",
    user_support_text=("trimmomatic is not bowtie"),
    citation_text=None,
)
plugin.methods.register_function(
    function=q2_trimmomatic.trim_paired,
    name="remove adpaters from paired sequences and quality trim",
    description="remove adpaters from paired sequences and quality trim",
    inputs={
        "paired_sequences": SampleData[PairedEndSequencesWithQuality],
    },
    parameters={"min_length": Int % Range(1, None)},
    outputs=[
        ("paired_end_trimmed", SampleData[PairedEndSequencesWithQuality]),
        ("unpaired_fwd", SampleData[SequencesWithQuality]),
        ("unpaired_rev", SampleData[SequencesWithQuality]),
    ],
    input_descriptions={"paired_sequences": "illumina paired sequence data"},
    parameter_descriptions={"min_length": "sequences shorter than this length will be discarded"},
    output_descriptions={
        "paired_end_trimmed": "paired trimmed sequence data",
        "unpaired_fwd": "unpaired trimmed sequence data fwd",
        "unpaired_rev": "unpaired trimmed sequence data rev",
    },
)
