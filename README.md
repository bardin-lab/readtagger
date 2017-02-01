tag_reads
_____

Tags reads in a BAM file based on other BAM files.

Installation
------

```
pip install https://github.com/bardin-lab/bamreadtagger/archive/master.zip
```

Useage
-------

To tag reads in file `a.bam` with file `b.bam` and output to path output.bam, type

```
tag_reads --tag_file a.bam --annotate_with b.bam ----output_file output.bam
```

This will by default tag reads with the RD, RR, MD and MR tags,
where the RD tag has detail mapping information for the current read,
while the MD tag has the information for the mate.
RR and MR contain the aligned reference (i.e chromosome).
The first letter can be changed on a per-file basis by appending ":first_letter_read:first_letter_mate" to the
file path. To change the above example into X for the read and Y for the mate, run:

```
tag_reads --tag_file a.bam --annotate_with b.bam:X:Z ----output_file output.bam
```

To tag one bam file using multiple alignment files, run:

```
tag_reads --tag_file a.bam --annotate_with b.bam:A:B c.bam:C:D ----output_file output.bam
```

Now reads that align in file `b.bam` will be tagged with AR, AD and BR, BD, while reads aligned in file `c.bam`
are marked with CR, CD and DR, DD.