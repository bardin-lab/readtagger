Readtagger
----------
.. image:: https://travis-ci.org/bardin-lab/readtagger.svg?branch=master
    :target: https://travis-ci.org/bardin-lab/readtagger

.. image:: https://coveralls.io/repos/github/bardin-lab/readtagger/badge.svg?branch=master
    :target: https://coveralls.io/github/bardin-lab/readtagger?branch=master

.. image:: https://badge.fury.io/py/readtagger.svg
    :target: https://badge.fury.io/py/readtagger

.. image:: https://anaconda.org/mvdbeek/readtagger/badges/version.svg
    :target: https://anaconda.org/mvdbeek/readtagger


Readtagger is a set of tools for describing the origin of reads in a sequenced genome. It can be used to verify and detect integration events of known, exogenous sequences (such as viruses) or endogenous sequences (such as transposons). In addition, there is the possibility of expanding this to de-novo detection without a priori knowledge of the inserted sequences).

Installation
------------

::

    pip install readtagger

Description
-----------

The readtagger module can be used as a tool for verifying transposable element (TE) insertions predicted by TEMP and TEPID, but it can also serve as the basis for detecting insertions using the findcluster command (see below). The readtagger module examines the origin of sequencing reads by aligning reads to multiple references and recording possible homologies as custom SAM tags. This works for single and paired-end reads, so that one obtains knowledge about a read that is aligned at a physical location in the genome as well as about its mate.


The readtagger command evaluates potential homologies, and (by default) discards homologies for alignments that do not improve the nucleotide usage, For example if the unclipped portion of an alignment has a homology with a transposable element, while we only care about the clipped, unaligned portion. Similarly, we discard homologies for aligned pairs that are well within the expected insert size. These precautions reduce noise in regions of the genome that have spurious homologies to transposable elements, as not all of these regions are masked by tools such as RepeatMasker.


To generate insertion calls, one can use the findcluster command.

findcluster operates in multiple phases, that are roughly:

1.	Divide the input file into approximately 100 KB bins. This allows parallelization over multiple CPU cores or machines. Splitting files may introduce artefacts when splitting in a cluster of tags, therefore we attempt to find regions that have no tagged alignments. 
2.	For each region we determine clusters of tagged alignments. A new cluster is created once additional alignments are not overlapping the current cluster coordinates. This is a rough first-pass that often includes multiple actual insertions and/or separates clusters that belong together, but does summarize all tagged alignments. Additionally we generate clusters of soft-clipped reads (independent of transposable element homology). 
3.	For each cluster, we ask whether there are any other clusters close by (the search area depends on the insert size distribution) that could potentially be merged because they provide evidence for the same insertion. All reachable clusters are tested by all-against-all alignments using edlib. Putative 5p and 3p sites are determined using a consensus alignment if split read evidence is available.
4.	We then ask whether the clusters are internally consistent. This means that alignments that point for example away from the putative insertion site are taken from the cluster and form a new, second cluster.
5.	We repeat steps 3 and 4 until the number of clusters at each pass stays constant.
6.	We assign soft-clipping clusters to transposable element clusters if any read in the soft clipping cluster is also member of a transposable element cluster. (This is important in a separate step that links transposable element insertion sites across samples).
7.	We then quantify the fragments that support a putative insertion (this includes alignments that did not reach the threshold for being tagged), the fragments that point against an insertion and the uninformative reads and we call genotypes and genotype likelihoods based on the evidence. Alignment evidence for an insertion and evidence against an insertion is tagged using separate tags, which allows for visual verification of the quantification. A VCF file describing the insertions, a BAM file containing relevant alignments and a fasta file of assembled sequences is written per 100 KB region. 
8.	Merge individual files


Usage
------

To tag reads in file ``a.bam`` with file ``b.bam`` and output to path
output.bam, type

::

    readtagger --tag_file a.bam --annotate_with b.bam ----output_file output.bam

This will by default tag reads with the AD, AR, BD and BR tags, where
the AD tag has detail mapping information for the current read, while
the BD tag has the information for the mate. AR and BR contain the
aligned reference (i.e chromosome). The first letter can be changed on a
per-file basis by appending ":first\_letter\_read:first\_letter\_mate"
to the file path. To change the above example into X for the read and Y
for the mate, run:

::

    readtagger --tag_file a.bam --annotate_with b.bam:X:Z ----output_file output.bam

To tag one bam file using multiple alignment files, run:

::

    readtagger --tag_file a.bam --annotate_with b.bam:A:B c.bam:C:D ----output_file output.bam

Now reads that align in file ``b.bam`` will be tagged with AR, AD and
BR, BD, while reads aligned in file ``c.bam`` are marked with CR, CD and
DR, DD.

Advanced usage
--------------

To see the advanced options, type:

::

    readtagger -h

Testing
-------

If you modify readtagger, you can run all tests by running tox:

::

    pip install tox
    tox

