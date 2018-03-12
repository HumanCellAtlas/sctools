import os
import tempfile
import pytest
import math
from sctools.metrics.gatherer import GatherGeneMetrics, GatherCellMetrics
import pandas as pd
import numpy as np

_data_dir = os.path.split(__file__)[0] + '/data'
_gene_sorted_bam = _data_dir + '/small-gene-sorted.bam'
_cell_sorted_bam = _data_dir + '/small-cell-sorted.bam'

_test_dir = tempfile.mkdtemp()
_gene_metric_output_file = _data_dir + '/gene_metrics.csv'
_cell_metric_output_file = _data_dir + '/cell_metrics.csv'
_testing_knowledge_scalar = pd.read_csv(
    _data_dir + '/small-gene-sorted_testing_knowledge_scalar.csv', index_col=0, squeeze=True,
    header=None)
_testing_knowledge_series = pd.read_csv(
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


gene_gatherer = GatherGeneMetrics(_gene_sorted_bam, _gene_metric_output_file)
gene_gatherer.extract_metrics()
_gene_metrics = pd.read_csv(_gene_metric_output_file, index_col=0)

cell_gatherer = GatherCellMetrics(_cell_sorted_bam, _cell_metric_output_file)
cell_gatherer.extract_metrics()
_cell_metrics = pd.read_csv(_cell_metric_output_file, index_col=0)


# def test_metrics_n_reads(gene_metrics):
#     reads_observed = gene_metrics['n_reads'].sum()
#     reads_expected = int(_testing_knowledge_scalar['n_reads'])
#     assert reads_expected == reads_observed


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 300),
    (_cell_metrics, 656),
])
def test_metrics_n_reads(metrics, expected_value):
    assert metrics['n_reads'].sum() == expected_value


# todo failing
def test_cell_metrics_mean_n_genes_observed():
    genes_observed = _cell_metrics['n_genes'].mean()
    assert math.isclose(genes_observed, 1.9827, abs_tol=1e-4), '%f != %f' % (genes_observed, 1.9827)


def test_gene_metrics_n_genes():
    genes_observed = _gene_metrics.shape[0]
    assert genes_observed == 8


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 88),
    (_cell_metrics, 249),
])
def test_metrics_n_molecules(metrics, expected_value):
    molecules_observed = metrics['n_molecules'].sum()
    assert molecules_observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 217),
    (_cell_metrics, 499),
])
def test_metrics_n_fragments(metrics, expected_value):
    fragments_observed = metrics['n_fragments'].sum()
    assert fragments_observed == expected_value


def test_metrics_highest_expression_gene():
    observed_max_gene = _gene_metrics['n_reads'].idxmax()
    expected_max_gene = _testing_knowledge_scalar['most_abundant']
    assert expected_max_gene == observed_max_gene


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 245),
    (_cell_metrics, 94),
])
def test_metrics_highest_read_count(metrics, expected_value):
    observed_max_gene_reads = metrics['n_reads'].max()
    assert observed_max_gene_reads == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 300),  # todo this is 100%, we should mangle a few in the testing data
    (_cell_metrics, 655),
])
def test_metrics_number_perfect_barcodes(metrics, expected_value):
    observed_perfect_barcodes = metrics['perfect_molecule_barcodes'].sum()
    assert observed_perfect_barcodes == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 300),  # todo this is 100%, should get some intronic or other reads
    (_cell_metrics, 609),
])
def test_reads_mapped_exonic(metrics, expected_value):
    observed = metrics['reads_mapped_exonic'].sum()
    assert observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 0),  # todo null case
    (_cell_metrics, 28),
])
def test_reads_mapped_intronic(metrics, expected_value):
    observed = metrics['reads_mapped_intronic'].sum()
    assert observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 0),  # todo null case
    (_cell_metrics, 19),
])
def test_reads_mapped_utr(metrics, expected_value):
    observed = metrics['reads_mapped_utr'].sum()
    assert observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 300),  # todo need to include at least 1 multi-mapper
    (_cell_metrics, 656),
])
def test_reads_mapped_uniquely(metrics, expected_value):
    observed = metrics['reads_mapped_uniquely'].sum()
    assert observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 90),
    (_cell_metrics, 107),
])
def test_duplicate_records(metrics, expected_value):
    observed = metrics['duplicate_reads'].sum()
    assert observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 29),
    (_cell_metrics, 2),
])
def test_spliced_reads(metrics, expected_value):
    observed = metrics['spliced_reads'].sum()
    assert observed == expected_value


# todo we need to further investigate why our assumption of duplicate_reads and fragments should
# sum to number of reads does not seem to be true in our test data
def test_relationship_of_duplicates_and_fragments():
    dup_and_fragments = _gene_metrics['duplicate_reads'].sum() + _gene_metrics['n_fragments'].sum()
    reads = _gene_metrics['n_reads'].sum()
    assert reads == dup_and_fragments


# def test_higher_order_metrics_by_gene():
#     higher_order_metrics = [
#         ('molecule_barcode_fraction_bases_above_30_mean'),
#         ('molecule_barcode_fraction_bases_above_30_variance'),
#         ('genomic_reads_fraction_bases_quality_above_30_mean'),
#         ('genomic_reads_fraction_bases_quality_above_30_variance'),
#         ('genomic_read_quality_mean'),
#         ('genomic_read_quality_variance'),
#         ('reads_per_molecule'),
#         ('reads_per_fragment'),
#         ('fragments_per_molecule')
#     ]
#     test_results = {}
#     for key in higher_order_metrics:
#         observed = np.nan_to_num(_gene_metrics[key].values)
#         expected = np.nan_to_num(_testing_knowledge_series[key].values)
#         test_results[key] = (np.allclose(observed, expected), observed, expected)
#     assert all(v[0] for v in test_results.values()), repr({k: v for k, v in test_results.items() if v[0] == False})


@pytest.mark.parametrize('metrics, key, expected_value', [
    (_cell_metrics, 'molecule_barcode_fraction_bases_above_30_mean',
     np.array([1.0000, 0.9500, 1.0000, 1.0000, 0.9778, 1.0000, 1.0000, 1.0000, 0.9833, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 0.9759, 1.0000, 1.0000, 0.9830, 1.0000, 1.0000,
               1.0000, 0.9778, 0.9783, 1.0000, 0.9800, 1.0000, 1.0000, 1.0000, 1.0000, 0.9500,
               1.0000, 0.9895, 1.0000, 0.9760, 1.0000, 1.0000, 1.0000, 0.9889, 1.0000, 0.9600,
               1.0000, 0.9909, 1.0000, 1.0000, 0.9556, 0.9800, 1.0000, 0.9000, 1.0000, 0.9588,
               1.0000, 1.0000, 0.9889, 0.8000, 0.9538, 0.9909, 0.9929, 0.9571])),
    (_cell_metrics, 'molecule_barcode_fraction_bases_above_30_variance',
     np.array(
         [np.nan, 0.0050, np.nan, np.nan, 0.0019, 0.0000, 0.0000, np.nan, 0.0015, np.nan, 0.0000,
          0.0000, np.nan, 0.0000, 0.0048, 0.0000, 0.0000, 0.0029, 0.0000, np.nan, 0.0000, 0.0044,
          0.0109, 0.0000, 0.0020, 0.0000, 0.0000, np.nan, 0.0000, 0.0100, np.nan, 0.0010, 0.0000,
          0.0052, 0.0000, 0.0000, 0.0000, 0.0011, 0.0000, 0.0162, 0.0000, 0.0016, 0.0000, np.nan,
          0.0178, 0.0020, np.nan, np.nan, 0.0000, 0.0163, np.nan, np.nan, 0.0011, np.nan, 0.0147,
          0.0018, 0.0007, 0.0306])),
    (_cell_metrics, 'genomic_reads_fraction_bases_quality_above_30_mean',
     np.array(
         [0.3980, 0.6786, 0.5000, 0.9796, 0.7800, 0.7811, 0.9337, 0.8469, 0.6743, 0.4565, 0.8622,
          0.9762, 0.4925, 0.7857, 0.7478, 0.8561, 0.6327, 0.7948, 0.8405, 0.4286, 0.7735, 0.6445,
          0.7291, 0.8520, 0.6711, 0.6123, 0.8238, 0.5000, 0.8376, 0.5137, 0.7526, 0.7584, 0.7574,
          0.8379, 0.8490, 0.5000, 0.5983, 0.7489, 0.7755, 0.8107, 0.6963, 0.8363, 0.8896, 0.6186,
          0.7549, 0.7151, 1.0000, 0.5306, 0.8347, 0.7340, 0.8367, 0.8878, 0.7347, 0.4592, 0.7718,
          0.7583, 0.8439, 0.7576])),
    (_cell_metrics, 'genomic_reads_fraction_bases_quality_above_30_variance',
     np.array(
         [np.nan, 0.1812, np.nan, np.nan, 0.0266, 0.0461, 0.0042, np.nan, 0.0387, np.nan, 0.0178,
          0.0000, np.nan, 0.0002, 0.0455, 0.0342, 0.0588, 0.0359, 0.0247, np.nan, 0.0400, 0.0436,
          0.0754, 0.0005, 0.1140, 0.0617, 0.0400, np.nan, 0.0230, 0.0491, np.nan, 0.0608, 0.0556,
          0.0367, 0.0215, 0.0860, 0.2182, 0.0564, 0.0008, 0.0395, 0.0330, 0.0433, 0.0063, np.nan,
          0.0366, 0.0778, np.nan, np.nan, 0.0114, 0.0391, np.nan, np.nan, 0.0193, np.nan, 0.0288,
          0.0444, 0.0311, 0.0558])),
    (_cell_metrics, 'genomic_read_quality_mean',
     np.array(
         [25.3776, 32.5051, 27.7755, 39.9184, 34.3639, 34.5969, 37.4592, 35.9490, 31.6345, 26.5870,
          36.7500, 39.5374, 28.0896, 33.7041, 33.6079, 36.2787, 30.8472, 34.8402, 35.9327, 24.7755,
          34.3603, 31.0934, 33.2880, 36.7092, 31.9647, 30.2158, 35.3956, 27.6837, 35.8674, 27.4527,
          34.3918, 33.7323, 33.6425, 35.9552, 35.5694, 27.4184, 30.0479, 33.4621, 34.6633, 35.2128,
          32.4619, 35.7690, 36.9963, 30.0722, 33.6353, 32.6708, 39.8721, 28.0510, 35.9388, 33.1278,
          35.8265, 36.6633, 32.7188, 26.6429, 34.1053, 34.0012, 36.0956, 33.7704])),
    (_cell_metrics, 'genomic_read_quality_variance',
     np.array([np.nan, 92.5078, np.nan, np.nan, 18.9818, 29.9521, 6.6724, np.nan, 25.4164, np.nan,
               12.8541, 0.3790, np.nan, 0.0019, 28.7815, 24.6669, 37.7402, 22.8765, 16.5399, np.nan,
               22.9679, 26.2414, 44.8249, 0.5740, 70.4607, 42.5318, 24.9536, np.nan, 14.0772,
               32.6389, np.nan, 38.1213, 34.4094, 23.2517, 13.9110, 48.9622, 117.2337, 32.9814,
               0.3850, 24.3135, 17.8765, 26.5847, 5.2099, np.nan, 22.5846, 48.2133, np.nan, np.nan,
               5.6775, 23.9395, np.nan, np.nan, 12.9322, np.nan, 18.1475, 29.6960, 20.7504,
               34.9055])),
    (_cell_metrics, 'reads_per_molecule',
     np.array([1.0000, 2.0000, np.nan, 1.0000, 9.0000, 2.4000, 2.0000, 1.0000, 3.0000, 1.0000,
               3.0000, 3.0000, 1.0000, np.nan, 2.4167, 4.3333, 1.2222, 5.8750, 1.3333, 1.0000,
               1.2000, 1.5000, 4.6000, 2.0000, 2.5000, 1.2000, 2.1429, 1.0000, 2.6364, 4.0000,
               1.0000, 2.1111, 1.7273, 6.2500, 5.0000, 1.3333, 2.0000, 2.2500, np.nan, 2.0000,
               4.3333, 3.9286, 2.2000, 1.0000, 1.5000, 1.6667, np.nan, 1.0000, 1.6667, 1.8889,
               1.0000, 1.0000, 2.2500, 1.0000, 9.7500, 11.0000, 4.0000, 1.5000])),
    (_cell_metrics, 'reads_per_fragment',
     np.array([1.0000, 1.0000, np.nan, 1.0000, 1.1250, 1.3333, 2.0000, 1.0000, 1.3333, 1.0000,
               1.2000, 3.0000, 1.0000, np.nan, 1.3182, 4.3333, 1.2222, 1.4688, 1.1429, 1.0000,
               1.2000, 1.2857, 1.5333, 2.0000, 1.2500, 1.2000, 1.2500, 1.0000, 1.3182, 1.0000,
               1.0000, 1.4615, 1.3571, 1.3889, 1.2500, 1.3333, 1.0000, 1.1250, np.nan, 1.1765,
               4.3333, 1.4474, 1.1000, 1.0000, 1.2857, 1.2500, np.nan, 1.0000, 1.2500, 1.3077,
               1.0000, 1.0000, 1.2857, 1.0000, 2.4375, 1.5714, 1.4737, 1.2353])),
    (_cell_metrics, 'fragments_per_molecule',
     np.array(
         [1.0000, 2.0000, np.nan, 1.0000, 8.0000, 1.8000, 1.0000, 1.0000, 2.2500, 1.0000, 2.5000,
          1.0000, 1.0000, np.nan, 1.8333, 1.0000, 1.0000, 4.0000, 1.1667, 1.0000, 1.0000, 1.1667,
          3.0000, 1.0000, 2.0000, 1.0000, 1.7143, 1.0000, 2.0000, 4.0000, 1.0000, 1.4444, 1.2727,
          4.5000, 4.0000, 1.0000, 2.0000, 2.0000, np.nan, 1.7000, 1.0000, 2.7143, 2.0000, 1.0000,
          1.1667, 1.3333, np.nan, 1.0000, 1.3333, 1.4444, 1.0000, 1.0000, 1.7500, 1.0000, 4.0000,
          7.0000, 2.7143, 1.2143])),
    (_gene_metrics, 'molecule_barcode_fraction_bases_above_30_mean',
     np.array([1.0000, 1.0000, 0.8000, 0.9885, 0.9833, 0.9857, 0.7000, 0.9444])),
    (_gene_metrics, 'molecule_barcode_fraction_bases_above_30_variance',
     np.array([np.nan, np.nan, np.nan, 0.0282, 0.0346, 0.0537, np.nan, 0.0849])),
    (_gene_metrics, 'genomic_reads_fraction_bases_quality_above_30_mean',
     np.array([0.8878, 0.3980, 0.4271, 0.8148, 0.7681, 0.7216, 0.1546, 0.5089])),
    (_gene_metrics, 'genomic_reads_fraction_bases_quality_above_30_variance',
     np.array([np.nan, np.nan, np.nan, 0.0282, 0.0346, 0.0537, np.nan, 0.0849])),
    (_gene_metrics, 'genomic_read_quality_mean',
     np.array([36.2143, 24.8469, 25.4792, 35.3664, 34.0956, 33.0364, 20.7423, 27.3078])),
    (_gene_metrics, 'genomic_read_quality_variance',
     np.array([np.nan, np.nan, np.nan, 18.4553, 21.6745, 33.6572, np.nan, 53.5457])),
    (_gene_metrics, 'reads_per_molecule',
     np.array([1.0000, 1.0000, 1.0000, 3.2500, 4.1525, 1.7500, 1.0000, 1.3846])),
    (_gene_metrics, 'reads_per_fragment',
     np.array([1.0000, 1.0000, 1.0000, 1.7333, 1.3920, 1.4000, 1.0000, 1.0588])),
    (_gene_metrics, 'fragments_per_molecule',
     np.array([1.0000, 1.0000, 1.0000, 1.8750, 2.9831, 1.2500, 1.0000, 1.3077])),
])
def test_higher_order_metrics_by_gene(metrics, key, expected_value):
    observed = sorted(np.nan_to_num(metrics[key].values).round(4))
    expected_value = sorted(np.nan_to_num(expected_value))
    assert observed == expected_value


@pytest.mark.parametrize('metrics, key, expected_value', [
    (_cell_metrics, 'n_molecules', 249),
    (_cell_metrics, 'n_fragments', 499),
    (_cell_metrics, 'fragments_with_single_read_evidence', 345),  # todo failing
    (_cell_metrics, 'molecules_with_single_read_evidence', 130),  # todo failing
    (_gene_metrics, 'n_molecules', 88),
    (_gene_metrics, 'n_fragments', 217),
    (_gene_metrics, 'fragments_with_single_read_evidence', 155),
    (_gene_metrics, 'molecules_with_single_read_evidence', 42),
])
def test_higher_order_metrics_across_genes(metrics, key, expected_value):
    observed = metrics[key].sum()
    assert observed == expected_value

