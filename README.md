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

This will by default tag reads with the AA, AR, MA and MR tags,
where the AA tag has detail mapping information for the current read,
while the MA tag has the information for the mate.
AR and MR contain the aigned reference (i.e chromosome).
The first letter of these tags canbe changed using the `--tag_prefix_self` and `--tag_prefix_mate` options.