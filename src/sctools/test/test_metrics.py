from io import StringIO
import os
import tempfile
import pytest
from sctools.metrics.gatherer import GatherGeneMetrics
import pandas as pd
import numpy as np

_data_dir = os.path.split(__file__)[0] + '/data'
_gene_sorted_bam = _data_dir + '/small-gene-sorted.bam'
_cell_sorted_bam = _data_dir + '/cell-sorted.bam'

_test_dir = tempfile.mkdtemp()
_gene_metric_output_file = _data_dir + '/gene_metrics.csv'
_cell_metric_output_file = _data_dir + '/cell_metrics.csv'
_testing_knowledge_scalar = pd.read_csv(
    _data_dir + '/small-gene-sorted_testing_knowledge_scalar.csv', index_col=0, squeeze=True,
    header=None)
_testing_knowledge_series = pd.DataFrame.from_csv(
    _data_dir + '/small-gene-sorted_testing_knowledge_series.csv')


# def test_calculate_cell_metrics():
#     platform.TenXV2.calculate_cell_metrics(
#         args=['-i', _cell_sorted_bam, '-o', _test_dir + '/gene_metrics.csv'])


# def test_calculate_gene_metrics():
#     platform.TenXV2.calculate_gene_metrics(
#         # args=['-i', _gene_sorted_bam, '-o', _test_dir + '/gene_metrics.csv'])
#         args=['-i', _gene_sorted_bam, '-o', 'gene_metrics.csv'])


# def test_merge_cell_metrics():
#     platform.TenXV2.merge_cell_metrics(
#         args=['-o', _test_dir + '/merged-cell-metrics.csv',
#               'cell_metrics_1.csv', 'cell_metrics_2.csv']
#     )


# def test_merge_gene_metrics():
#     platform.TenXV2.merge_gene_metrics(
#         args=['-o', _test_dir + '/merged-gene-metrics.csv',
#               'gene_metrics_1.csv', 'gene_metrics_2.csv']
#     )


@pytest.fixture(scope='module')
def gene_metrics():
    """run gene metrics on the testing file and pull in the outputs"""
    gene_gatherer = GatherGeneMetrics(_gene_sorted_bam, _gene_metric_output_file)
    gene_gatherer.extract_metrics()
    return pd.read_csv(_gene_metric_output_file, index_col=0)


def test_metrics_file_identifies_correct_n_reads(gene_metrics):
    reads_observed = gene_metrics['n_reads'].sum()
    reads_expected = int(_testing_knowledge_scalar['n_reads'])
    assert reads_expected == reads_observed


def test_metrics_file_identifies_correct_n_genes(gene_metrics):
    genes_observed = gene_metrics.shape[0]
    genes_expected = int(_testing_knowledge_scalar['n_genes'])
    assert genes_expected == genes_observed


def test_metrics_file_identifies_correct_n_molecules(gene_metrics):
    molecules_observed = gene_metrics['n_molecules'].sum()
    molecules_expected = int(_testing_knowledge_scalar['n_molecules'])
    assert molecules_expected == molecules_observed


def test_metrics_file_identifies_correct_n_fragments(gene_metrics):
    fragments_observed = gene_metrics['n_fragments'].sum()
    fragments_expected = int(_testing_knowledge_scalar['n_fragments'])
    assert fragments_expected == fragments_observed


def test_metrics_file_identifies_correct_highest_expression_gene(gene_metrics):
    observed_max_gene = gene_metrics['n_reads'].idxmax()
    expected_max_gene = _testing_knowledge_scalar['most_abundant']
    assert expected_max_gene == observed_max_gene


def test_metrics_file_identifies_correct_highest_expression_count(gene_metrics):
    observed_max_gene_reads = gene_metrics['n_reads'].max()
    expected_max_gene_reads = int(_testing_knowledge_scalar['most_abundant_gene_n_observations'])
    assert expected_max_gene_reads == observed_max_gene_reads


def test_metrics_file_identifies_correct_number_perfect_barcodes(gene_metrics):
    observed_perfect_barcodes = gene_metrics['perfect_molecule_barcodes'].sum()
    expected_perfect_barcodes = int(_testing_knowledge_scalar['perfect_molecule_barcodes'])
    assert observed_perfect_barcodes == expected_perfect_barcodes


def test_reads_mapped_exonic(gene_metrics):
    observed = gene_metrics['reads_mapped_exonic'].sum()
    expected = int(_testing_knowledge_scalar['reads_mapped_exonic'])
    assert observed == expected


# todo this is a zero-case in the current test data; need some intronic reads
def test_reads_mapped_intronic(gene_metrics):
    observed = gene_metrics['reads_mapped_intronic'].sum()
    expected = int(_testing_knowledge_scalar['reads_mapped_intronic'])
    assert observed == expected


# todo this is a zero-case in the current test data; need some utr reads
def test_reads_mapped_utr(gene_metrics):
    observed = gene_metrics['reads_mapped_utr'].sum()
    expected = int(_testing_knowledge_scalar['reads_mapped_utr'])
    assert observed == expected


def test_reads_mapped_uniquely(gene_metrics):
    observed = gene_metrics['reads_mapped_uniquely'].sum()
    expected = int(_testing_knowledge_scalar['reads_mapped_uniquely'])
    assert observed == expected


def test_duplicate_records(gene_metrics):
    observed = gene_metrics['duplicate_reads'].sum()
    expected = int(_testing_knowledge_scalar['duplicate_reads'])
    assert observed == expected


    # todo this test currently fails because we believe our data is not deduped correctly, need to run on data run on latest optimus pipeline
def test_relationship_of_duplicates_and_fragments(gene_metrics):
    dup_and_fragments = gene_metrics['duplicate_reads'].sum() + gene_metrics['n_fragments'].sum()
    reads = gene_metrics['n_reads'].sum()
    assert reads == dup_and_fragments


def test_spliced_reads(gene_metrics):
    observed = gene_metrics['spliced_reads'].sum()
    expected = int(_testing_knowledge_scalar['spliced_reads'])
    assert observed == expected


def test_higher_order_metrics_by_gene(gene_metrics):
    higher_order_metrics = [
        ('molecule_barcode_fraction_bases_above_30_mean'),
        ('molecule_barcode_fraction_bases_above_30_variance'),
        ('genomic_reads_fraction_bases_quality_above_30_mean'),  # These four tests are commented out because
        ('genomic_reads_fraction_bases_quality_above_30_variance'), # the notebook equivalent dont ignore soft clipped bases
        ('genomic_read_quality_mean'),
        ('genomic_read_quality_variance'),
        ('reads_per_molecule'),
        ('reads_per_fragment'),
        ('fragments_per_molecule')
    ]
    test_results = {}
    for key in higher_order_metrics:
        observed = np.nan_to_num(gene_metrics[key].values)
        expected = np.nan_to_num(_testing_knowledge_series[key].values)
        test_results[key] = (np.allclose(observed, expected), observed, expected)
    assert all(v[0] for v in test_results.values()), repr({k: v for k, v in test_results.items() if v[0] == False})


def test_higher_order_metrics_across_genes(gene_metrics):
    higher_order_metrics = [
        ('n_molecules', int),
        ('n_fragments', int),
        ('fragments_with_single_read_evidence', int),
        ('molecules_with_single_read_evidence', int)
    ]
    test_results = {}
    for (key, func) in higher_order_metrics:
        observed = gene_metrics[key].sum()
        expected = func(_testing_knowledge_scalar[key])
        test_results[key] = (observed == expected, observed, expected)
    assert all(v[0] for v in test_results.values()), repr({k: v for k, v in test_results.items() if v[0] == False})

