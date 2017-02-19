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

Tags reads in a BAM file based on other BAM files.

Installation
------------

::

    pip install readtagger

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

