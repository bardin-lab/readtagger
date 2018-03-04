import os
import tempfile

import six


def write_sequences(sequences, output_path=None, tmp_dir=None):
    """
    Take sequences and write them out to output_path.

    >>> import shutil
    >>> tmp_dir = tempfile.mkdtemp(prefix='readtagger_write_sequences')
    >>> d = {'A': 'ATGC'}
    >>> l = ['ATGC', 'ATGC']
    >>> s = 'ATGC'
    >>> paths = [write_sequences(sequences=x) for x in [d, l, s]]
    >>> assert len(paths) == 3
    >>> shutil.rmtree(tmp_dir, ignore_errors=True)
    """
    fd = None
    if not output_path:
        fd, output_path = tempfile.mkstemp(dir=tmp_dir)
    with open(output_path, 'w') as out:
        if isinstance(sequences, dict):
            for qname, seq in sequences.items():
                out.write(">%s\n%s\n" % (qname, seq))
        elif isinstance(sequences, list):
            for qname, seq in enumerate(sequences):
                out.write(">%s\n%s\n" % (qname, seq))
        elif isinstance(sequences, six.string_types):
            out.write(">0\n%s\n" % sequences)
    if fd:
        os.close(fd)
    return output_path
