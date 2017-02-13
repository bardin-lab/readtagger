import subprocess
import pysam


class BamAlignmentWriter(object):
    """
    Simple wrapper around pysam.AlignmentFile that uses sambamba for multithreaded compressed writing
    """
    def __init__(self, path, template=None, header=None, threads=2):
        self.template = template
        self.header= header
        self.args = ['sambamba', 'view', '-S', '-f', 'bam', '-t', "%s" % threads, '/dev/stdin' , '-o', path]
        self.proc = subprocess.Popen(self.args, stdin=subprocess.PIPE)

    def write(self, r):
        return self.af.write(r)

    def close(self):
        self.af.close()
        self.proc.stdin.close()

    def __enter__(self):
        self.af = pysam.AlignmentFile(self.proc.stdin, mode="wh", template=self.template, header=self.header)
        return self.af

    def __exit__(self, type, value, traceback):
        self.close()


class BamAlignmentReader(object):
    """
    Simple wrapper around pysam.AlignmentFile that uses sambamba for reading if input file is a bam file
    """
    def __init__(self, path):
        self.bam = pysam.AlignmentFile(path).is_bam
        self.path = path
        self.args = ['sambamba', 'view', '-h', path]

    def close(self):
        self.af.close()
        if self.bam:
            self.proc.stdout.close()

    def __enter__(self):
        if self.bam:
            self.proc = subprocess.Popen(self.args, stdout=subprocess.PIPE)
            self.af = pysam.AlignmentFile(self.proc.stdout)
        else:
            self.af = pysam.AlignmentFile(self.path)
        return self.af

    def __exit__(self, type, value, traceback):
        self.close()