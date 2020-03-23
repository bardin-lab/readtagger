from setuptools import setup

__VERSION__ = '0.5.22'

ENTRY_POINTS = '''
        [console_scripts]
        add_matesequence=readtagger.cli.add_matesequence:annotate_mate
        allow_dovetailing=readtagger.cli.allow_dovetailing:allow_dovetailing
        annotate_softclipped_reads=readtagger.cli.annotate_softclipped_reads:annotate_softclipped_reads
        confirm_insertions=readtagger.cli.classify_somatic_insertions:confirm_insertions
        findcluster=readtagger.cli.findcluster:findcluster
        merge_clusterfinder_vcfs=readtagger.cli.merge_findcluster_vcf:merge_findcluster
        plot_coverage=readtagger.cli.plot_coverage:plot_coverage
        pysamtools_view=readtagger.cli.pysamtools_view_cli:pysamtools_view
        readtagger=readtagger.cli.readtagger_cli:readtagger
        update_mapq=readtagger.cli.update_mapq:update_mapq
        write_supplementary_fastq=readtagger.cli.write_supplementary_fastq:write_supplementary_fastq
        extract_variants=readtagger.cli.extract_variants:extract_variants
        normalize_readsizes=readtagger.cli.normalize_readsizes:normalize_readsizes
        summarize_fragments=readtagger.cli.summarize_fragments:cli
'''

requirements = ['bcbio-gff',
                'biopython',
                'cached_property',
                'click',
                'contextlib2',
                'compare-reads',
                'edlib',
                'matplotlib',  # for plotting coverage, refactor into separate package
                'mappy',
                'multiprocessing_logging',
                'pandas',
                'pysam',
                'scipy',
                'six']

readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

setup(
    name='readtagger',
    version=__VERSION__,
    packages=['readtagger', 'readtagger.cli'],
    description="Tags reads in a BAM file based on other BAM files.",
    long_description=readme + '\n\n' + history,
    install_requires=requirements,
    entry_points=ENTRY_POINTS,
    keywords='Bioinformatics',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Environment :: Console',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    extras_require={
        'testing': ["pytest", "pytest-datadir", "tox", "planemo", "cookiecutter", "bumpversion"],
    },
    url='https://github.com/bardin-lab/readtagger',
    license='MIT',
    author='Marius van den Beek',
    author_email='m.vandenbeek@gmail.com',
)
