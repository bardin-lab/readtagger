import os
import tempfile

import six


def write_sequences(sequences, output_path=None, tmp_dir=None):
    """Take sequences and write them out to output_path."""
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
