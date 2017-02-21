.. :changelog:

History
-------

.. to_doc

---------------------
0.2.0(2017-02-21)
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

