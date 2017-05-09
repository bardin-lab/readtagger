import click

from readtagger.readtagger import TagManager
from readtagger import VERSION


def parse_file_tags(filetags):
    """
    Parse list of filetags from commandline.

    :param filetags: list of strings with filepath.
                     optionally appended by the first letter that should be used for read and mate
    :return: annotate_with, tag_prefix, tag_prefix_mate

    >>> filetags = ('file_a:A:B', 'file_b:C:D', 'file_c')
    >>> annotate_with, tag_prefix, tag_prefix_mate = parse_file_tags(filetags)
    >>> annotate_with == ['file_a', 'file_b', 'file_c'] and tag_prefix == ['A', 'C', 'A'] and tag_prefix_mate == ['B', 'D', 'B']
    True
    >>>
    """
    annotate_with = []
    tag_prefix = []
    tag_prefix_mate = []
    if not isinstance(filetags, tuple):
        filetags = (filetags,)
    for filetag in filetags:
        if ':' in filetag:
            filepath, tag, tag_mate = filetag.split(':')
            annotate_with.append(filepath)
            tag_prefix.append(tag.upper())
            tag_prefix_mate.append(tag_mate.upper())
        else:
            annotate_with.append(filetag)
            tag_prefix.append('A')  # Default is A for read, B for mate
            tag_prefix_mate.append('B')
    return annotate_with, tag_prefix, tag_prefix_mate


@click.command()
@click.option('-t', '--target_path', help="Annotate reads in this file.", required=True, type=click.Path(exists=True))
@click.option('-s', '--source_path',
              help="Tag reads in target_path if reads are aligned in these files."
                   "Append `:A:B` to tag first letter of tag describing read as A, "
                   "and first letter of tag describing the mate as B",
              multiple=True,
              required=True)
@click.option('-o', '--output_path', help="Write bam file to this path", required=True)
@click.option('-r', '--reference_fasta', help='Reference fasta to align clipped reads to', default=None, required=False)
@click.option('--allow_dovetailing/--no_allow_dovetailing', default=True,
              help="Sets the proper pair flag (0x0002) to true if reads dovetail [reads reach into or surpass the mate sequence].")
@click.option('--discard_if_proper_pair/--no_discard_if_proper_pair', default=True, help="Discard an alternative flag if the current read is in a proper pair.")
@click.option('--discard_suboptimal_alternate_tags/--no_discard_suboptimal', default=True,
              help="By default cigarstrings of the alternative tags are compared and alternates that are not explaining the current cigar "
                   "strings are discarded. Use this option to keep the alternative tags (effectively restoring the behaviour of readtagger < 0.1.4)")
@click.option('--discarded_path', help="Write reads with discarded tags to this path", type=click.Path(exists=False), required=False)
@click.option('--verified_path', help="Write reads with verified tags to this path", type=click.Path(exists=False), required=False)
@click.option('--cores', default=1, help='Number of cores to use for tagging reads.')
@click.version_option(version=VERSION)
def main(**kwargs):
    """Tag reads in alignment file `target_path` with reads in `source_path`."""
    source_path, tag_prefix_self, tag_prefix_mate = list(zip(*parse_file_tags(kwargs.pop('source_path'))))[0]
    kwargs['source_path'] = source_path
    kwargs['tag_prefix_self'] = tag_prefix_self
    kwargs['tag_prefix_mate'] = tag_prefix_mate
    TagManager(**kwargs)
