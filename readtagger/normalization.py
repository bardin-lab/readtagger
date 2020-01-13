"""Normalize length of multiple long read fastq files."""
import gzip
import bz2

from collections import (
    defaultdict,
    namedtuple,
)
from itertools import zip_longest
from operator import (
    attrgetter,
    neg,
)

from Bio import SeqIO
from sortedcontainers import SortedKeyList

# if after splitting less than REMAINING_LENGTH_LIMIT nuncleotides are left discard that fragment
REMAINING_LENGTH_LIMIT = 500
IndexToReadlength = namedtuple('IndexToReadlength', 'index readlength')


def open_by_suffix(filename, mode='rt'):
    """Open compressed or uncompressed files."""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    elif filename.endswith('.bz2'):
        return bz2.BZ2file(filename, mode)
    else:
        return open(filename, mode[0])


def calculate_split(lists, result=None):
    """
    Given multiple lists with read lengths truncate reads so that reads are equally long.

    Returns a list of lists. The outer list is is a new list of truncated readlengths,
    the inner list contains IndexToReadlength instances in the order of ``lists``.
    """
    result = result or []
    aux_lists = defaultdict(lambda: SortedKeyList(key=lambda x: neg(x.readlength)))
    for list_index in range(len(lists)):
        # initialize aux lists
        aux_lists[list_index].add(IndexToReadlength(index=-1, readlength=-1))
    for t in zip_longest(*lists, fillvalue=IndexToReadlength(index=-1, readlength=-1)):
        # when iterator is exhausted index_to_readlength instance will have index and readlength at -1
        # t is a tuple with items at the current iteration
        lengths = []
        for list_index, index_to_readlength in enumerate(t):
            if aux_lists[list_index]:
                aux_length = aux_lists[list_index].pop()
                if aux_length.readlength > index_to_readlength.readlength:
                    if index_to_readlength.readlength > 0:
                        aux_lists[list_index].add(index_to_readlength)
                    lengths.append(aux_length)
                else:
                    if aux_length.readlength > 0:
                        aux_lists[list_index].add(aux_length)
                    if index_to_readlength.readlength > 0:
                        lengths.append(index_to_readlength)
            else:
                lengths.append(index_to_readlength)
        min_length = min(lengths, key=attrgetter('readlength'))
        intermediate_result = []
        for list_index, old_length in enumerate(lengths):
            remaining_length = old_length.readlength - min_length.readlength
            if remaining_length > REMAINING_LENGTH_LIMIT:
                aux_lists[list_index].add(IndexToReadlength(index=old_length.index, readlength=remaining_length))
            intermediate_result.append(IndexToReadlength(index=old_length.index, readlength=min_length.readlength))
        result.append(intermediate_result)
    return result


def flatten_calculated_split_lists(split_result):
    """Flatten nested lists of readlength tuples per file."""
    r = defaultdict(list)
    for inner_list in split_result:
        for i, readlength in enumerate(inner_list):
            r[i].append(readlength)
    for k, v in r.items():
        v.sort(key=attrgetter('index'), reverse=True)  # so pop gives us the first index
    return r


def get_read_sizes(path):
    """Create list of IndexToReadlength tuples sorted by readlength."""
    length_index = []
    with open_by_suffix(path) as fh:
        for i, (_, sequence, _) in enumerate(SeqIO.QualityIO.FastqGeneralIterator(fh)):
            length_index.append(IndexToReadlength(index=i, readlength=len(sequence)))
    return sorted(length_index, key=attrgetter('readlength'), reverse=True)


def split_reads(path, final_sizes, output_path):
    """Given a list of target sizes and read indexes split reads in ``path`` and write split reads out to ``output_path``."""
    final_size = final_sizes.pop()
    while final_size.index < 0:
        final_size = final_sizes.pop()
    with open_by_suffix(path) as fh_in, open_by_suffix(output_path, 'wt') as fh_out:
        for i, (title, sequence, quality) in enumerate(SeqIO.QualityIO.FastqGeneralIterator(fh_in)):
            try:
                while final_size.index == i:
                    record = "@%s\n%s\n+\n%s\n" % (title, sequence[:final_size.readlength], quality[:final_size.readlength])
                    fh_out.write(record)
                    sequence = sequence[final_size.readlength:]
                    quality = quality[final_size.readlength:]
                    final_size = final_sizes.pop()
                if final_size.index < i:
                    # get next final size
                    final_size = final_sizes.pop()
                elif final_size.index > i:
                    # get next record
                    continue
            except IndexError as e:
                if 'pop from empty list' not in str(e):
                    raise


def split_fastq_files(input_paths, output_paths):
    """
    "Normalize length of multiple long read fastq files.

    Establishes a list of read size tuples for each input fastq, where each list is sorted from longest to shortest read.
    The shorter length of the longest read length across all lists will then be the target read length.
    Fragments of the longer reads will be stored in an adjacent list so that these may be re-used when comparing shorter
    longest reads.
    """
    read_sizes = [get_read_sizes(f) for f in input_paths]
    r = flatten_calculated_split_lists(calculate_split(read_sizes))
    for input_path, final_size, output_path in zip(input_paths, r.values(), output_paths):
        split_reads(path=input_path, final_sizes=final_size, output_path=output_path)
