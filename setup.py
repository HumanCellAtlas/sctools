from setuptools import setup

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Natural Language :: English",
    "License :: OSI Approved :: BSD License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: Implementation :: PyPy",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name='sctools',
    version='0.1.6',
    description='Utilities for large-scale distributed single cell data processing',
    url='https://github.com/humancellatlas/sctools.git',
    author='Ambrose J. Carr',
    author_email='mail@ambrosejcarr.com',
    package_dir={'': 'src'},
    packages=['sctools', 'sctools/test', 'sctools/metrics'],
    install_requires=[
        'pysam',
        'numpy',
        'pandas',
        'pytest',
    ],
    entry_points={
            'console_scripts': [
                'Attach10xBarcodes = sctools.platform:TenXV2.attach_barcodes',
                'SplitBam = sctools.platform:GenericPlatform.split_bam',
                'CalculateGeneMetrics = sctools.platform:GenericPlatform.calculate_gene_metrics',
                'CalculateCellMetrics = sctools.platform:GenericPlatform.calculate_cell_metrics',
                'MergeGeneMetrics = sctools.platform:GenericPlatform.merge_gene_metrics',
                'MergeCellMetrics = sctools.platform:GenericPlatform.merge_cell_metrics',
            ]
    },
    classifiers=CLASSIFIERS,
    include_package_data=True,
)
