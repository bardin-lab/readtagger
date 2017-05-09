from .helpers import reference_fasta  # noqa: F401
from readtagger.tag_softclip import TagSoftClip

COMPLEX = 'extended_annotated_updated_all_reads.bam'


def test_tag_softclip(datadir, tmpdir, reference_fasta):  # noqa: D103, F811
    input_path = datadir[COMPLEX]
    output_bam = tmpdir.join('output.bam').strpath
    TagSoftClip(source=input_path, output_path=output_bam, reference_fasta=reference_fasta)
