from setuptools import setup

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Natural Language :: English",
    "License :: OSI Approved :: BSD License",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="sctools",
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    description="Utilities for large-scale distributed single cell" +
                "data processing",
    url="https://github.com/humancellatlas/sctools.git",
    author="Ambrose J. Carr",
    author_email="mail@ambrosejcarr.com",
    package_dir={"": "src"},
    packages=["sctools", "sctools/test", "sctools/metrics"],
    install_requires=[
        "gffutils",
        "numpy",
        "pandas",
        "pysam",
        "pytest",
        "pytest-cov",
        "sphinx",
        "sphinxcontrib-websupport",
        "sphinx_rtd_theme",
        "setuptools_scm>=3.1.0",
        "setuptools>=40.4.3",
        "scipy>=1.0.0",
        "crimson>=0.3.0",
    ],
    entry_points={
        "console_scripts": [
            "AttachBarcodes = sctools.platform:BarcodePlatform." +
            "attach_barcodes",
            "Attach10xBarcodes = sctools.platform:TenXV2.attach_barcodes",
            "SplitBam = sctools.platform:GenericPlatform.split_bam",
            "CalculateGeneMetrics = sctools.platform:GenericPlatform." +
            "calculate_gene_metrics",
            "CalculateCellMetrics = sctools.platform:GenericPlatform." +
            "calculate_cell_metrics",
            "MergeGeneMetrics = sctools.platform:GenericPlatform." +
            "merge_gene_metrics",
            "MergeCellMetrics = sctools.platform:GenericPlatform." +
            "merge_cell_metrics",
            "CreateCountMatrix = sctools.platform:GenericPlatform." +
            "bam_to_count_matrix",
            "MergeCountMatrices = sctools.platform:GenericPlatform." +
            "merge_count_matrices",
            "TagSortBam = sctools.platform:GenericPlatform.tag_sort_bam",
            "VerifyBamSort = sctools.platform:GenericPlatform.verify_bam_sort",
            "GroupQCs = sctools.platform:GenericPlatform.group_qc_outputs",
        ]
    },
    classifiers=CLASSIFIERS,
    include_package_data=True,
)
