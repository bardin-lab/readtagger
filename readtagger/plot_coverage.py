import itertools
import os
from collections import (
    defaultdict,
    Mapping
)
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402


# from https://gist.github.com/angstwad/bf22d1822c38a92ec0a9
def dict_merge(dct, merge_dct):
    """
    Merge dicts recursively.

    Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    """
    for k in merge_dct:
        if (k in dct and isinstance(dct[k], dict) and isinstance(merge_dct[k], Mapping)):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]


def dd():
    """Return defaultdict with dict factory."""
    return defaultdict(dict)


def get_coverage(file, label, regions=None, nth=1, readcount=-1):
    """Get coverage for every `nth` position from alignment file."""
    readcount = float(readcount)
    contigs_coverage = defaultdict(dd)
    with pysam.AlignmentFile(file) as f:
        if isinstance(regions, str):
            regions = [regions]
        for region in regions:
            for pileup_pos in f.pileup(region=region, max_depth=20000):
                if pileup_pos.pos % nth == 0:
                    reference_name = f.get_reference_name(pileup_pos.tid)  # Seems to have broken in pysam 0.1.14
                    contigs_coverage[reference_name][label][pileup_pos.reference_pos] = pileup_pos.nsegments / (readcount / 10**6) if readcount else 0
    return contigs_coverage


def plot_coverage(contigs_coverage, style='ggplot', nrows=8):
    """Plot coverage in contigs_coverage."""
    plt.style.use(style)
    figs = []
    i = 0
    for title, data in contigs_coverage.items():
        if i % nrows == 0:
            fig_axes = plt.subplots(nrows=nrows, figsize=(10, 20))
            fig = fig_axes[0]
            axes = fig_axes[1:]
            plt.subplots_adjust(hspace=0.3)
            figs.append(fig)
        i += 1
        nrow = i - (int(i / nrows) * nrows)
        ax = pd.DataFrame.from_dict(data).plot(kind='area',
                                               title=title,
                                               stacked=False,
                                               alpha=0.5,
                                               fig=fig,
                                               ax=axes[0][nrow])
        ax.legend(bbox_to_anchor=(1.1, 1), loc="upper right")
    return figs


def mp_get_coverage(args):
    """Wrap get_coverage for multiprocessing Pool implementation."""
    return get_coverage(*args)


def get_total_coverage(file, cores=1):
    """Get number of reads in file."""
    return int(pysam.view('-c', "-@%d" % cores, file).strip())


def plot_coverage_in_regions(files, labels, output_path, regions=None, cores=1, total_reads=None):
    """
    Plot coverage for `files`, where files are multiple BAM files.

    `output_path` is the path at which the plot should be saved
    `regions` can be speficied, these should be chromosome names for now.
    """
    if not regions:
        regions = pysam.AlignmentFile(files[0]).references
    for f in files:
        if not os.path.exists("%s.bai" % f):
            pysam.index(f)
    if not total_reads:
        total_reads = [get_total_coverage(file, cores=cores) for file in files]
    starmap_args = [(file, label, region, 10, reads) for (file, label, reads), region in itertools.product(zip(files, labels, total_reads), regions)]
    if cores == 1:
        r = itertools.starmap(get_coverage, starmap_args)
    else:
        pool = ProcessPoolExecutor(max_workers=cores)
        r = pool.map(mp_get_coverage, starmap_args)
        pool.shutdown()
    contigs_coverage = next(r)
    for d in r:
        dict_merge(contigs_coverage, d)
    figs = plot_coverage(contigs_coverage)
    if output_path:
        with PdfPages(output_path) as pdf:
            for f in figs:
                pdf.savefig(f)
    return figs
