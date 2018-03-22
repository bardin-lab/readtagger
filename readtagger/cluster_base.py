import logging
import os
from itertools import chain
from cached_property import cached_property

from .gff_io import (
    sort_gff,
    write_gff_cluster
)

logger = logging.getLogger(__name__)


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
            if hasattr(self, 'softclip_finder'):
                clusters = chain(self.clusters, self.softclip_finder.clusters)
            else:
                clusters = self.clusters
            write_cluster(clusters=clusters,
                          header=self.header,
                          output_path=self.output_gff,
                          sample=self.sample_name,
                          threads=self.threads)
            if self.threads < 2:
                sort_gff(self.output_gff, output_path=self.output_gff)
