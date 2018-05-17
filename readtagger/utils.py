import os
import shlex
import subprocess
import tempfile
from .instance_lru import lru_cache

COMPLEMENTARY_SEQUENCES = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


@lru_cache(maxsize=10000)
def revcom(string):
    """Build reverse complement of string."""
    return "".join([COMPLEMENTARY_SEQUENCES[s] for s in string[::-1]])


def sort(input_path, output_path, sort_cmd="sort -k 1,1 -k4,4n"):
    """Sort text file at input path."""
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(output_path), prefix="sort_")
    header_lines = []
    with open(input_path) as gff_in, open(tmp, 'w') as out:
        for line in gff_in:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                out.write(line)
        out.close()
    with open(output_path, 'w') as out:
        cmd = shlex.split(sort_cmd)
        cmd.append(tmp)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, close_fds=True)
        # No need to add newlines, those should be present already
        out.write("".join(header_lines))
        for line in p.stdout:
            out.write(line.decode())
    os.close(fd)


def overlap(start1, end1, start2, end2, tolerance=0):
    """Check that range (start1, end1) overlaps with (start2, end2)."""
    # Taken from https://nedbatchelder.com/blog/201310/range_overlap_in_two_compares.html
    return end1 + tolerance >= start2 and end2 + tolerance >= start1
