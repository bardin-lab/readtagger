.. :changelog:

History
-------

.. to_doc


---------------------
0.3.8 (2017-03-14)
---------------------
* Fix multithreaded blast runs and do not fail on empty Blast output

  ---------------------
0.3.7 (2017-03-14)
---------------------
* Add reference_fasta for blast functionality in findcluster galaxy wrapper
* Cleanup BlastProcessor object at the end of a run
* New test fixture for downloading reference fasta
* Implement first draft of blasting cluster contigs
* Limit clusterfinding to reads with `AD` or `BD` tag

---------------------
0.3.6 (2017-03-09)
---------------------
* Use query_alignment_sequence instead of query_sequence for compatibility with BWA MEMs -Y option
* Keep suppl. reads without primary reads by default

---------------------
0.3.5 (2017-03-08)
---------------------
* Install blast in travis environemnt and use conda-compatible virtualenv
* Add tools to extract supplementary alignments and tool to update MAPQ scores of supplementary alignments
* Add blast module
* Move to click as commandline argument parser for find cluster and write_supplementary_fastq
* Refactor code organization in project.
* Update ROADMAP.rst

---------------------
0.3.4 (2017-03-02)
---------------------
* Speed up assembly steps using multithreading
* Implement a cache for the Cluster.can_join method

---------------------
0.3.3 (2017-03-02)
---------------------
* Fix a crash when writing GFF for a cluster of hardclipped reads
* Change confusing variable names and copypasted docstring

---------------------
0.3.2 (2017-03-02)
---------------------
* Fix another crash when tuple starts with 1,2,7 or 8

---------------------
0.3.1 (2017-03-02)
---------------------
* Fix a crash when a mismatch is the last item in a cigartuple

---------------------
0.3.0 (2017-03-02)
---------------------

* Add a galaxy tool for the findcluster script
* Add new script that finds clusters of reads and outputs GFF or BAM files with these clusters.
* Implement writing clusters as GFF files
* Implement writing out reads with cluster number annotated in CD tag.
* Implement merging of clusters based on whether reads contribute to common contigs
* Use cached-property where it makes sense
* Add module to find, join and annotate clusters of reads
* Represent cigartuple as namedtuple
* Add a Roadmap file
* Add more logic for finding ends of insertions and
* Manipulate cluster of reads to find TSDs
* Add module for cap3 assembly and manipulation of assembled reads
* Fix conda recipe script entrypoints

---------------------
0.2.0 (2017-02-21)
---------------------
* Reformat help text in galaxy wrappers
* Add add_matesequence script to add the sequence of the mate of the current read as a tag
* Add option to discard alternative tag if read is a proper pair
* Stitch cigars that are separated by I or D events
* Add a tag tuple that knows how to format itself
* Update README.rst example with current default tag prefix
* Test with and without discarding verified reads
* Symlink test-files that are shared with the galaxy test, add testcase for allow_dovetailing script
* Fix HISTORY.rst formatting

---------------------
0.1.13(2017-02-17)
---------------------
* Add instructions for development
* Install planemo in deployment step

---------------------
0.1.12(2017-02-17)
---------------------
* Test deployment again

---------------------
0.1.11 (2017-02-17)
---------------------
* Test deployment

---------------------
0.1.10 (2017-02-17)
---------------------
* Fix toolshed deployment

---------------------
0.1.9 (2017-02-17)
---------------------
* Add automated deployment to Galaxy Toolshed
* Add instructions for development and release process

---------------------
0.1.8 (2017-02-17)
---------------------
* Minor release to test release process

---------------------
0.1.7 (2017-02-17)
---------------------
* Extend testing with coverage testing
* Automate deployment to pypi and conda
* Register project with pyup.io

---------------------
0.1.6 (2017-02-16)
---------------------
* Rename to readtagger
* Fix bug with stdin closing file descriptor too early, leading to corrupt
  BAM files
* Extend testing

---------------------
0.1.5 (2017-02-12)
---------------------
* Add option (-wd) to write suboptimal tag into separate BAM file
* Add option (-wv) to write verified tags into separate BAM file
* Performance improvments by letting sambamba handle BAM reading
  and writing. Also elimininate regualr expression to parse cigarstring

---------------------
0.1.4 (2017-02-10)
---------------------
* Add option (-k) to keep alternative tags if they do not
  explain the softclipped read any better.
  Default is to discard them.

---------------------
0.1.3.2 (2017-02-08)
---------------------
* Fix dovetailing script

---------------------
0.1.3 (2017-02-07)
---------------------
* Add option to allow dovetailing in alignment files when tagging reads
* Add separate entrypoint for standalone script

---------------------
0.1.2 (2017-02-05)
---------------------
* Add conda recipe
* Python3 string fix

---------------------
0.1.0 (2017-02-05)
---------------------
* Initial version
