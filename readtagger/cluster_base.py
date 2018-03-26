import logging
import os
from itertools import chain
from cached_property import cached_property

from .gff_io import (
    sort_gff,
    write_gff_cluster
)
from .vcf_io import write_vcf

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


class ToOutput(object):
    """Provides logic for writing clusters and softclip clusters."""

    def to_output(self, output_path, write_func, sort_func):
        """Write clusters as GFF file."""
        logger.info("Writing clusters of GFF (%s)", self.region or 0)
        if output_path:
            if hasattr(self, 'softclip_finder'):
                clusters = chain(self.clusters, self.softclip_finder.clusters)
            else:
                clusters = self.clusters
            write_func(clusters=clusters,
                       header=self.header,
                       output_path=output_path,
                       sample_name=self.sample_name,
                       threads=self.threads)
            if self.threads < 2:
                sort_func(input_path=output_path, output_path=output_path)


class ToGffMixin(ToOutput):
    """Provide a `to_gff` function."""

    def to_gff(self, output_path):
        """Write clusters as GFF file."""
        self.to_output(output_path, write_func=write_gff_cluster, sort_func=sort_gff)


class ToVcfMixin(ToOutput):
    """Provide a `to_vcf` function."""

    def to_vcf(self, output_path):
        """Write clusters as VCF file."""
        self.to_output(output_path=output_path, write_func=write_vcf, sort_func=sort_gff)
