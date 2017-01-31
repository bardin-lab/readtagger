try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

requirements = ['pysam','six']

ENTRY_POINTS = '''
        [console_scripts]
        tag_reads=tag_reads.tag_reads:main
'''

setup(
    name='tag_reads',
    version='0.1.0',
    packages=['tag_reads'],
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
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
    ],
    url='https://github.com/bardin-lab/tag_reads',
    license='MIT',
    author='Marius van den Beek',
    author_email='m.vandenbeek@gmail.com',
    description='Tags reads in BAM files based on alignments in additional BAM files.'
)
