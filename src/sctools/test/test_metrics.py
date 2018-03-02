import os
import tempfile
from itertools import product

import pytest
import pandas as pd
import numpy as np

from sctools.metrics.gatherer import GatherGeneMetrics, GatherCellMetrics

_data_dir = os.path.split(__file__)[0] + '/data'
_gene_sorted_bam = _data_dir + '/small-gene-sorted.bam'
_cell_sorted_bam = _data_dir + '/small-cell-sorted.bam'

_test_dir = tempfile.mkdtemp()
_gene_metric_output_file = _data_dir + '/gene_metrics.csv'
_cell_metric_output_file = _data_dir + '/cell_metrics.csv'
_scalar_gene_testing_knowledge = pd.read_csv(
    _data_dir + '/small-gene-sorted_testing_knowledge_scalar.csv', index_col=0, squeeze=True,
    header=None)
_scalar_cell_testing_knowledge = pd.read_csv(
    _data_dir + '/small-cell-sorted_testing_knowledge_scalar.csv', index_col=0, squeeze=True,
    header=None)


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


cell_gatherer = GatherCellMetrics(_cell_sorted_bam, _cell_metric_output_file)
cell_gatherer.extract_metrics()
_cell_metrics = pd.read_csv(_cell_metric_output_file, index_col=0)

gene_gatherer = GatherGeneMetrics(_gene_sorted_bam, _gene_metric_output_file)
gene_gatherer.extract_metrics()
_gene_metrics = pd.read_csv(_gene_metric_output_file, index_col=0)

parameters = [
    # ['n_genes', int, np.sum],  # data['genes_detected'] for cells
    ['n_molecules', int, np.sum],
    ['n_fragments', int, np.sum],
    # ('most_abundant', str, lambda x: x.index[np.argmax(x, 0]]],
    # ['most_abundant_gene_n_observations', int, np.max],
    ['perfect_molecule_barcodes', int, np.sum],
    ['reads_mapped_exonic', int, np.sum],
    ['reads_mapped_intronic', int, np.sum],  # todo test data has zero of these
    ['reads_mapped_utr', int, np.sum],  # todo test data has zero of these
    ['reads_mapped_uniquely', int, np.sum],
    ['duplicate_reads', int, np.sum],
    ['spliced_reads', int, np.sum],
]

metrics = [
    [_cell_metrics, _scalar_cell_testing_knowledge],
    [_gene_metrics, _scalar_gene_testing_knowledge],
]

test_combinations = [m + p for m, p in product(metrics, parameters)]


@pytest.mark.parametrize("observed,expected", [
    (_cell_metrics, _scalar_cell_testing_knowledge),
    (_gene_metrics, _scalar_gene_testing_knowledge),
])
def test_n_reads(observed, expected):
    reads_observed = observed['n_reads'].sum()
    reads_expected = int(expected['n_reads'])
    assert reads_expected == reads_observed


@pytest.mark.parametrize("observed,expected,index,converter,reduce_func", test_combinations)
def test_scalar_metrics(observed, expected, index, converter, reduce_func):
    obs_val = reduce_func(observed[index])
    exp_val = converter(expected[index])
    assert obs_val == exp_val, 'observed: %f != expected: %f for %s' % (obs_val, exp_val, index)


# gene-specific metric tests

def test_n_genes_genes():
    genes_observed = _gene_metrics.shape[0]
    genes_expected = int(_scalar_gene_testing_knowledge['n_genes'])
    assert genes_expected == genes_observed


# cell-specific metrics
def test_n_genes_cells():
    genes_observed = np.sum(_cell_metrics['genes_detected'])
    genes_expected = int(_scalar_cell_testing_knowledge['genes_detected'])
    assert genes_expected == genes_observed


# def test_highest_expression_gene(gene_metrics):
#     observed_max_gene = gene_metrics['n_reads'].idxmax()
#     expected_max_gene = _testing_knowledge['most_abundant']
#     assert expected_max_gene == observed_max_gene
#
#
# def test_highest_expression_count(gene_metrics):
#     observed_max_gene_reads = gene_metrics['n_reads'].max()
#     expected_max_gene_reads = int(_testing_knowledge['most_abundant_gene_n_observations'])
#     assert expected_max_gene_reads == observed_max_gene_reads


# def test_relationship_of_duplicates_and_fragments(gene_metrics):
#     dup_and_fragments = gene_metrics['duplicate_reads'].sum() + gene_metrics['n_fragments'].sum()
#     reads = gene_metrics['n_reads'].sum()
#     assert reads == dup_and_fragments

@pytest.mark.parametrize('index,converter', [
    ('molecule_barcode_fraction_bases_above_30_mean', float),
    ('molecule_barcode_fraction_bases_above_30_variance', float),
    ('genomic_reads_fraction_bases_quality_above_30_mean', float),
    ('genomic_reads_fraction_bases_quality_above_30_variance', float),
    ('genomic_read_quality_mean', float),
    ('genomic_read_quality_variance', float),
    ('n_molecules', int),
    ('n_fragments', int),
    ('reads_per_molecule', float),
    ('reads_per_fragment', float),
    ('fragments_per_molecule', float),
    ('fragments_with_single_read_evidence', int),
    ('molecules_with_single_read_evidence', int),
])
def test_higher_order_metrics(index, converter):
    observed = _gene_metrics[index].mean()
    expected = converter(_scalar_gene_testing_knowledge[index])
    assert observed == expected, '%s failed' % index
