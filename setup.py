from setuptools import setup

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Natural Language :: English",
    "License :: OSI Approved :: BSD License",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: Implementation :: PyPy",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name='sctools',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    description='Utilities for large-scale distributed single cell data processing',
    url='https://github.com/humancellatlas/sctools.git',
    author='Ambrose J. Carr',
    author_email='mail@ambrosejcarr.com',
    package_dir={'': 'src'},
    packages=['sctools', 'sctools/test', 'sctools/metrics'],
    install_requires=[
        'gffutils',
        'numpy',
        'pandas',
        'pysam',
        'pytest',
        'pytest-cov',
        'sphinx',
        'sphinxcontrib-napoleon',
        'sphinx_rtd_theme',
        'setuptools_scm'
        'scipy>=1.0.0',
    ],
    entry_points={
            'console_scripts': [
                'Attach10xBarcodes = sctools.platform:TenXV2.attach_barcodes',
                'SplitBam = sctools.platform:GenericPlatform.split_bam',
                'CalculateGeneMetrics = sctools.platform:GenericPlatform.calculate_gene_metrics',
                'CalculateCellMetrics = sctools.platform:GenericPlatform.calculate_cell_metrics',
                'MergeGeneMetrics = sctools.platform:GenericPlatform.merge_gene_metrics',
                'MergeCellMetrics = sctools.platform:GenericPlatform.merge_cell_metrics',
                'CreateCountMatrix = sctools.platform:GenericPlatform.bam_to_count_matrix',
                'MergeCountMatrices = sctools.platform:GenericPlatform.merge_count_matrices',
            ]
    },
    classifiers=CLASSIFIERS,
    include_package_data=True,
)
