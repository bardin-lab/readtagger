import subprocess
import pysam


class BamAlignmentWriter(object):
    """
    Simple wrapper around pysam.AlignmentFile that uses sambamba for multithreaded compressed writing
    """
    def __init__(self, outpath, template, threads=4):
        self.template = template
        self.args = ['sambamba', 'view', '-S', '-f', 'bam', '-t', "%s" % threads, '/dev/stdin' , '-o', outpath]
        self.proc = subprocess.Popen(self.args, stdin=subprocess.PIPE)

    def write(self, r):
        return self.af.write(r)

    def close(self):
        self.af.close()
        self.proc.stdin.close()

    def __enter__(self):
        self.af = pysam.AlignmentFile(self.proc.stdin, mode="wh", template=self.template)
        return self.af

    def __exit__(self, type, value, traceback):
        self.close()
