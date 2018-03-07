import logging
import os
from cached_property import cached_property

from .gff_io import (
    write_cluster,
    sort_gff
)


class SampleNameMixin(object):
    """Provide a sample name property."""

    @cached_property
    def sample_name(self):
        """Return sample name if passed in manually, else guess sample name from input file."""
        if not self._sample_name:
            basename = os.path.basename(self.input_path)
            if '.' in basename:
                basename = basename.rsplit('.', 1)[0]
            return basename
        else:
            return self._sample_name


class ToGffMixin(object):
    """Provide a `to_gff` function."""

    def to_gff(self):
        """Write clusters as GFF file."""
        logging.info("Writing clusters of GFF (%s)", self.region or 0)
        if self.output_gff:
            write_cluster(clusters=self.clusters,
                          header=self.header,
                          output_path=self.output_gff,
                          sample=self.sample_name,
                          threads=self.threads)
            if self.threads < 2:
                sort_gff(self.output_gff, output_path=self.output_gff)
