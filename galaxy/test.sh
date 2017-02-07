#!/usr/bin/env bash

set -e

if grep -v 'python tag_reads.py' bam_tag_reads.xml
then
    sed -i.bak 's/tag_reads -t/python \$__tool_directory__\/tag_reads.py -t/g' bam_tag_reads.xml
fi
cp ../tag_reads/tag_reads.py .
planemo test --conda_prefix ~/miniconda3
