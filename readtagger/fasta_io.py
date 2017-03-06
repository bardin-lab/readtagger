def write_sequences(sequences, output_path):
    """Take sequences and write them out to output_path."""
    with open(output_path, 'w') as out:
        for qname, seq in sequences.items():
            out.write(">%s\n%s\n" % (qname, seq))
