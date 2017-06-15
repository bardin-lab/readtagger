import tempfile

import six


def write_sequences(sequences, output_path=None, tmp_dir=None):
    """Take sequences and write them out to output_path."""
    if not output_path:
        _, output_path = tempfile.mkstemp(dir=tmp_dir)
    with open(output_path, 'w') as out:
        if isinstance(sequences, dict):
            for qname, seq in sequences.items():
                out.write(">%s\n%s\n" % (qname, seq))
        elif isinstance(sequences, list):
            for qname, seq in enumerate(sequences):
                out.write(">%s\n%s\n" % (qname, seq))
        elif isinstance(sequences, six.string_types):
            out.write(">0\n%s\n" % sequences)
    return output_path
