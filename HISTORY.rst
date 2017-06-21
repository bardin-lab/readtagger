.. :changelog:

History
-------

.. to_doc

---------------------
0.3.25 (2017-06-21)
---------------------
* Refine cluster coordinates using an Assembly strategy
* Fix GFF sorting on python 3
* Improve BWA alignment settings (default to intractg plus -Y) and add align_contigs method to SimpleAligner
* Add pysamtools_view command
* Improve cluster-splitting
* Add multiprocessing-logging recipe
* Only output BWA stderr if the exit code is not zero
* Add a function to sort gff files
* Close open file descriptors
* Make imprecise insertion sites more realistic
* Fix read_index property
* Adapt readtagger to higher coverage datasets
* Fix readtagger crash when not producing discard tag file.
* Add number of mates for left and right support to GFF
* Split clusters that start with reverse reads conatining only BD tags

---------------------
0.3.24 (2017-05-11)
---------------------
* Split cluster if there are multiple polarity switches between Forward and Reverse orientation
* Manipulate copy of cigarlist to avoid numpy issue

---------------------
0.3.23 (2017-05-09)
---------------------
* Expose reference fasta option in bam_readtagger.xml

---------------------
0.3.22 (2017-05-09)
---------------------
* Move readtagger CLI form argparse to click
* Index bamfile if neccesary
* Replace multipocessing pool with ProcessPoolExecutor
* Set the matesequence while tagging reads
* Fix false positives in readtagger module
* Do cap3 assembly in shared memory if passing --shm_dir or if SHM_DIR environment variable is defined
* Parallelize findlcluster by splitting input bam
* Add check_call.py script for rapidly verifying IGV screenshots

---------------------
0.3.21 (2017-04-27)
---------------------
* Fix crash when determining reference name

---------------------
0.3.20 (2017-04-27)
---------------------
* Guess the best TE match and write it into GFF Parent
* Fix case where input files are already sorted
* Remove blast from requirements

---------------------
0.3.19 (2017-04-27)
---------------------
* Skip creating tempdirs in current working directory
* Remove blast-specific files
* Switch to using BWA for annotating detected insertions
* Add more logging and default to not changing sort order unless specifically demanded
* Do dovetailing on coordinate-sorted file

---------------------
0.3.18 (2017-04-25)
---------------------
* Fix small outputs due to switching of `-t` and `-a` options

---------------------
0.3.17 (2017-04-25)
---------------------
* Fix file seeking
* Update dependencies

---------------------
0.3.16 (2017-04-23)
---------------------
* Parallelize readtagger

---------------------
0.3.15 (2017-04-20)
---------------------
* Do not count reads as support if both AD and BD tag contribute to an insertion
* Remove sambamba support

---------------------
0.3.14 (2017-04-19)
---------------------
* Perform readtagging on readname sorted files.
* Catch possible errors
* Add BWA alignment module to replace Blast

---------------------
0.3.13 (2017-04-05)
---------------------
* Add possibility to output cluster contigs as fasta

---------------------
0.3.12 (2017-03-31)
---------------------
* Fix and accelerate the calculation of nref (=non support evidence)
* Update priors and genotype frequrencies to a more realistic model

---------------------
0.3.11 (2017-03-28)
---------------------
* Add a testcase for genotyping module
* Stream over full alignment file instead of fetching regions,
  pysam.AlignmentFile.fetch is too slow

---------------------
0.3.10 (2017-03-26)
---------------------
* Revert local conda dependency resolution
* Fix readtagger.add_mate to work also if one mate is unmapped

---------------------
0.3.9 (2017-03-26)
---------------------
* Add a genotyping module
* Keep tags for alternative alignments if mates are not in a proper pair

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
