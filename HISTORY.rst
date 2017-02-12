.. :changelog:

History
-------

.. to_doc

---------------------
0.1.5 (2017-02-12)
---------------------
* Add option (-wd) to write suboptimal tag into separate BAM file
* Add option (-wv) to write verified tags into separate BAM file
* Performance improvemtns by letting sambamba handle BAM reading
  and writing. Also elimininate regualr expression to parse cigarstring.

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

