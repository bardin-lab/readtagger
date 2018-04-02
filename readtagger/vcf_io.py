import datetime
import os
from string import Template

import pysam
import pysam.bcftools

from .utils import sort

try:
    from readtagger import VERSION
except Exception:
    VERSION = '1.2.3'
try:
    from tempfile import TemporaryDirectory
except ImportError:
    from backports.tempfile import TemporaryDirectory

# It is possible to define sample to genome mappings as shown below:
# #META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
# #META=<ID=Disease,Type=String,Number=.,Values=[None, Cancer]>
# #META=<ID=Ethnicity,Type=String,Number=.,Values=[AFR, CEU, ASN, MEX]>
# #META=<ID=Tissue,Type=String,Number=.,Values=[Blood, Breast, Colon, Lung, ?]>
# #SAMPLE=<ID=Sample1,Assay=WholeGenome,Ethnicity=AFR,Disease=None,Description="Patient germline genome from unaffected",DOI=url>
# #SAMPLE=<ID=Sample2,Assay=Exome,Ethnicity=CEU,Disease=Cancer,Tissue=Breast,Description="European patient exome from breast cancer">


# # 1.4.9 Pedigree field format
# # It is possible to record relationships between genomes using the following syntax:
# ##PEDIGREE=<ID=TumourSample,Original=GermlineID>
# ##PEDIGREE=<ID=SomaticNonTumour,Original=GermlineID>
# ##PEDIGREE=<ID=ChildID,Father=FatherID,Mother=MotherID>
# ##PEDIGREE=<ID=SampleID,Name_1=Ancestor_1,...,Name_N=Ancestor_N>

CONTIG_TEMPLATE = Template("""##contig=<ID=$contig_name,length=$contig_length>""")
DATE = datetime.datetime.now()
SOURCE = "readtagger-v{version}".format(version=VERSION)
FILE_DATE = "{year}{month}{day}".format(year=DATE.year, month=DATE.month, day=DATE.day)
VCF_HEADER_LINE = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'

VCF_HEADER_TEMPLATE = Template("""##fileformat=VCFv4.1
##fileDate=$date
##source=$source
$contigs
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Maximum MAPQ of evidence supporting the variant">
##INFO=<ID=VALID_TSD,Number=0,Type=Flag,Description="Insertion is flanked by a Target Site Duplication.">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DEL:ME,Description="Deletion of mobile element present in reference">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=INS:ME,Description="Insertion of mobile element">
##ALT=<ID=SOFTCLIP,Description="Softclipped alignments">
##ALT=<ID=SOFTCLIP:5P,Description="Softclipped alignments where softclipping removes the 5p side">
##ALT=<ID=SOFTCLIP:3P,Description="Softclipped alignments where softclipping removes the 3p side">
##FORMAT=<ID=MENAME,Number=1,Type=String,Description="Mobile element name">
##FORMAT=<ID=MESTART,Number=1,Type=Integer,Description="Mobile element reference start">
##FORMAT=<ID=MEEND,Number=1,Type=Integer,Description="Mobile element reference end">
##FORMAT=<ID=MEASSEMBLY5,Number=.,Type=String,Description="Sequences assembled at 5p breakpoint of mobile element insertion">
##FORMAT=<ID=MEASSEMBLY3,Number=.,Type=String,Description="Sequences assembled at 3p breakpoint of mobile element insertion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">
##FORMAT=<ID=SU5,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant from the 5p side">
##FORMAT=<ID=SU3,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant form the 3p side">
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant">
##FORMAT=<ID=PE5,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant from the 5p side">
##FORMAT=<ID=PE3,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant from the 3p side">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of split reads supporting the variant">
##FORMAT=<ID=SR5,Number=1,Type=Integer,Description="Number of split reads supporting the variant from the 5p side">
##FORMAT=<ID=SR3,Number=1,Type=Integer,Description="Number of split reads supporting the variant from the 3p side">
##FORMAT=<ID=MSP,Number=1,Type=Integer,Description="Number of fragments supporting insertion with both mates. Indicates short insertions.">
##FORMAT=<ID=CLIP_CONSENSUS,Number=1,Type=String,Description="Consensus sequence for SOFTCLIP event. The consensus is listed in the 5p to 3p direction relative to the alignment.">
$header_line
""")  # noqa: E501


def get_vcf_contig_lines(header):
    """Return properly formatted contig lines for VCF/BCF headers."""
    template = CONTIG_TEMPLATE
    contigs = (template.substitute(contig_name=n, contig_length=l) for n, l in zip(header.references, header.lengths))
    return "\n".join(contigs)


def get_vcf_header(header, sample_name):
    """Return complete VCF header."""
    contigs = get_vcf_contig_lines(header)
    return VCF_HEADER_TEMPLATE.substitute(date=DATE,
                                          source=SOURCE,
                                          contigs=contigs,
                                          header_line="%s\t%s\n" % (VCF_HEADER_LINE, sample_name))


def write_vcf(output_path, clusters, header, sample_name, **kwargs):
    """Write clusters as VCF."""
    header_content = get_vcf_header(header=header, sample_name=sample_name)
    with TemporaryDirectory(prefix="tmp_vcf_header") as temp_dir:
        vcf_header_tmp = os.path.join(temp_dir, 'header.vcf')
        with open(vcf_header_tmp, 'w') as header_out:
            header_out.write(header_content)
        with pysam.VariantFile(vcf_header_tmp) as vcf_sample_file:
            vcf_header = vcf_sample_file.header
        with pysam.VariantFile(output_path, 'w', header=vcf_header) as vcf_out:
            for cluster in clusters:
                record = vcf_out.new_record()
                for k, v in cluster.vcf_mandatory.items():
                    v = getattr(cluster, v)
                    setattr(record, k, v)
                for k, v in cluster.vcf_info.items():
                    record.info[k] = getattr(cluster, v)
                for k, v in cluster.vcf_sample.items():
                    if isinstance(v, list):
                        v = [getattr(cluster, _) for _ in v]
                    else:
                        v = getattr(cluster, v)
                    record.samples[sample_name][k] = v
                vcf_out.write(record)


def merge_vcf_files(vcf_files, output_path, sort_output=True):
    """Merge vcf files."""
    # Ideally we'd be able to use
    # pysam.bcftools.merge('-o', output_path, *vcf_files)
    # But this only works for bgzipped vcf files.
    all_vfs = [pysam.VariantFile(vcf_file) for vcf_file in vcf_files if os.path.exists(vcf_file)]
    try:
        header = all_vfs[0].header
        for vf in all_vfs:
            header.merge(vf.header)
        with pysam.VariantFile(output_path, 'w', header=header) as output:
            for vf in all_vfs:
                for r in vf:
                    output.write(r)
    finally:
        for vf in all_vfs:
            vf.close()
    if sort_output:
        sort_vcf(input_path=output_path, output_path=output_path)


def sort_vcf(input_path, output_path):
    """Sort VCF."""
    return sort(input_path=input_path, output_path=output_path, sort_cmd="sort -k 1,1 -k2,2n")
