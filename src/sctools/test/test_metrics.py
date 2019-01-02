import fileinput
import math
import os
import tempfile
from typing import Callable

import numpy as np
import pandas as pd
import pytest
from sctools.metrics.gatherer import GatherGeneMetrics, GatherCellMetrics, MetricGatherer
from sctools.metrics.merge import MergeCellMetrics, MergeGeneMetrics
from sctools.platform import TenXV2

"""
Testing Data Definition & Acquisition

The data hardcoded into this file come from the two notebooks associated with these metrics:
characterize-cell-testing-data.ipynb and characterize-gene-testing-data.ipynb. In these
notebooks, the testing .bam files are loaded into memory and interrogated for each of the
metrics in question using pandas and numpy commands. These independent calculation provide the 
hard-coded data found in these tests. When testing data is changed, the notebook can be updated 
to re-calculate the values found in this file.
"""

# set the input and output directories, using a tempdir to automatically clean up generated files
_data_dir = os.path.split(__file__)[0] + '/data'
_test_dir = tempfile.mkdtemp()
os.makedirs(_test_dir, exist_ok=True)

# note, to inspect these testing files, please install samtools and use the following command:
# samtools view <filename> | less

# set the input files
_gene_sorted_bam = os.path.join(_data_dir, 'small-gene-sorted.bam')
_cell_sorted_bam = os.path.join(_data_dir, 'small-cell-sorted.bam')
_cell_sorted_bam_missing_cell_barcodes = os.path.join(_data_dir, 'cell-sorted-missing-cb.bam')

# specify filenames for temporary metrics outputs that are used in the following tests
_gene_metric_output_file = os.path.join(_test_dir, 'gene_metrics.csv.gz')
_cell_metric_output_file = os.path.join(_test_dir, 'cell_metrics.csv.gz')
_cell_metric_output_file_missing_cell_barcodes = os.path.join(_test_dir, 'cell_metrics_missing_cb.csv.gz')

# run the gene metrics suite
gene_gatherer = GatherGeneMetrics(_gene_sorted_bam, _gene_metric_output_file)
gene_gatherer.extract_metrics()
_gene_metrics = pd.read_csv(_gene_metric_output_file, index_col=0)

# run the cell metrics suite
cell_gatherer = GatherCellMetrics(_cell_sorted_bam, _cell_metric_output_file)
cell_gatherer.extract_metrics()
_cell_metrics = pd.read_csv(_cell_metric_output_file, index_col=0)

# run the cell metrics suite
cell_gatherer_missing_cbs = GatherCellMetrics(_cell_sorted_bam_missing_cell_barcodes, _cell_metric_output_file_missing_cell_barcodes)
cell_gatherer_missing_cbs.extract_metrics()
_cell_metrics_missing_cbs = pd.read_csv(_cell_metric_output_file_missing_cell_barcodes, index_col=0)


def test_calculate_cell_metrics_cli():
    """test the sctools cell metrics CLI invocation"""
    cell_metrics_csv = os.path.join(_test_dir, 'cell_metrics.csv')
    return_call = TenXV2.calculate_cell_metrics(
        args=['-i', _cell_sorted_bam, '-o', cell_metrics_csv])
    assert return_call == 0


def test_calculate_gene_metrics_cli():
    """test the sctools gene metrics CLI invocation"""
    gene_metrics_csv = os.path.join(_test_dir, 'gene_metrics.csv')
    return_call = TenXV2.calculate_gene_metrics(
        args=['-i', _gene_sorted_bam, '-o', gene_metrics_csv])
    assert return_call == 0


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 300),
    (_cell_metrics, 656),
])
def test_metrics_n_reads(metrics, expected_value):
    """test that the metrics identify the correct read number"""
    assert metrics['n_reads'].sum() == expected_value


def test_cell_metrics_mean_n_genes_observed():
    """
    test that the GatherCellMetrics method identifies the correct number of genes per cell, on
    average.
    """
    genes_observed = _cell_metrics['n_genes'].mean()
    assert math.isclose(genes_observed, 1.9827, abs_tol=1e-4), '%f != %f' % (genes_observed, 1.9827)


def test_gene_metrics_n_genes():
    """Test that GatherGeneMetrics identifies the total number of genes in the test file"""
    genes_observed = _gene_metrics.shape[0]
    assert genes_observed == 8


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 88),
    (_cell_metrics, 249),
])
def test_metrics_n_molecules(metrics, expected_value):
    """Test that each metric identifies the total number of molecules in the test file

    Molecules are defined as a unique combination of {cell barcode, molecule barcode, gene}
    """
    molecules_observed = metrics['n_molecules'].sum()
    assert molecules_observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 217),
    (_cell_metrics, 499),
])
def test_metrics_n_fragments(metrics, expected_value):
    """Test that each metric identifies the total number of fragments in the test file.

    Fragments are defined as a unique combination of {cell barcode, molecule barcode, strand,
    position, chromosome}
    """
    fragments_observed = metrics['n_fragments'].sum()
    assert fragments_observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 'AL627309.7'),
    (_cell_metrics, 'AAACCTGGTAGAAGGA'),
])
def test_metrics_highest_expression_class(metrics, expected_value):
    """
    for gene metrics, this is the highest expression gene. For cell metrics, this is the highest
    expression cell.
    """
    observed_max_gene = metrics['n_reads'].idxmax()
    assert observed_max_gene == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 245),
    (_cell_metrics, 94),
])
def test_metrics_highest_read_count(metrics, expected_value):
    """
    Test that each metric identifies the what the highest read count associated with any single
    entity
    """
    observed_max_gene_reads = metrics['n_reads'].max()
    assert observed_max_gene_reads == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 300),  # todo this is 100%, we should mangle a few in the testing data
    (_cell_metrics, 655),
])
def test_metrics_number_perfect_molecule_barcodes(metrics, expected_value):
    """Test that each metric correctly identifies the number of perfect molecule barcodes where UB == UR"""
    observed_perfect_barcodes = metrics['perfect_molecule_barcodes'].sum()
    assert observed_perfect_barcodes == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_cell_metrics, 650),
    (_cell_metrics_missing_cbs, 12861)
])
def test_metrics_number_perfect_cell_barcodes(metrics, expected_value):
    """Test that each metric correctly identifies the number of perfect cell barcodes where CB == CR"""
    observed_perfect_cell_barcodes = metrics['perfect_cell_barcodes'].sum()
    assert observed_perfect_cell_barcodes == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 300),  # todo this is 100%, should get some intronic or other reads
    (_cell_metrics, 609),
])
def test_reads_mapped_exonic(metrics, expected_value):
    """Test that each metric identifies the number of reads mapped to an exon (XF=='CODING')"""
    observed = metrics['reads_mapped_exonic'].sum()
    assert observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 0),  # todo null case
    (_cell_metrics, 28),
])
def test_reads_mapped_intronic(metrics, expected_value):
    """Test that each metric identifies the number of reads mapped to an intron (XF=='INTRONIC')"""
    observed = metrics['reads_mapped_intronic'].sum()
    assert observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 0),  # todo null case
    (_cell_metrics, 19),
])
def test_reads_mapped_utr(metrics, expected_value):
    """Test that each metric identifies the number of reads mapped to a UTR (XF=='UTR')"""
    observed = metrics['reads_mapped_utr'].sum()
    assert observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 300),  # todo need to include at least 1 multi-mapper
    (_cell_metrics, 656),
])
def test_reads_mapped_uniquely(metrics, expected_value):
    """Uniquely mapping reads will be tagged with NH==1"""
    observed = metrics['reads_mapped_uniquely'].sum()
    assert observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 90),
    (_cell_metrics, 107),
])
def test_duplicate_records(metrics, expected_value):
    """Duplicate records are identified by the 1024 bit being set in the sam flag"""
    observed = metrics['duplicate_reads'].sum()
    assert observed == expected_value


@pytest.mark.parametrize('metrics, expected_value', [
    (_gene_metrics, 29),
    (_cell_metrics, 2),
])
def test_spliced_reads(metrics, expected_value):
    """
    This pipeline defines spliced reads as containing an N segment of any length in the cigar string
    """
    observed = metrics['spliced_reads'].sum()
    assert observed == expected_value


# todo failing
# @pytest.mark.parametrize('metrics', [_gene_metrics, _cell_metrics])
# def test_relationship_of_duplicates_and_fragments(metrics):
#     """
#     We expect the number of duplicates and fragments to add up to the total number of reads. The
#     rationale is that any read that is not a duplicate should be a distinct fragment, under our
#     definitions.
#
#     This fails because of (1) N-base and 2-base cell barcode correction errors and (2)
#     fragment calculationes currently do not account for soft clipping. Fixing these will cause
#     this test to pass
#     """
#     dup_and_fragments = metrics['duplicate_reads'].sum() + metrics['n_fragments'].sum()
#     reads = metrics['n_reads'].sum()
#     assert reads == dup_and_fragments


@pytest.mark.parametrize('metrics', [_gene_metrics, _cell_metrics])
def test_fragments_number_is_greater_than_molecule_number(metrics):
    """
    There should always be more fragments than molecules, as the minimum definition of a molecule is
    a fragment covered by a single read
    """
    assert np.all(metrics['n_molecules'] >= 1)
    assert np.all(metrics['n_fragments'] >= 1)
    assert np.all(metrics['n_fragments'] >= metrics['n_molecules'])


@pytest.mark.parametrize('metrics, key, expected_value', [
    (_cell_metrics, 'molecule_barcode_fraction_bases_above_30_mean',
     np.array([1.0000, 0.9500, 1.0000, 1.0000, 0.9778, 1.0000, 1.0000, 1.0000, 0.9833, 1.0000,
               1.0000, 1.0000, 1.0000, 1.0000, 0.9759, 1.0000, 1.0000, 0.9830, 1.0000, 1.0000,
               1.0000, 0.9778, 0.9783, 1.0000, 0.9800, 1.0000, 1.0000, 1.0000, 1.0000, 0.9500,
               1.0000, 0.9895, 1.0000, 0.9760, 1.0000, 1.0000, 1.0000, 0.9889, 1.0000, 0.9600,
               1.0000, 0.9909, 1.0000, 1.0000, 0.9556, 0.9800, 1.0000, 0.9000, 1.0000, 0.9588,
               1.0000, 1.0000, 0.9889, 0.8000, 0.9538, 0.9909, 0.9929, 0.9571])),
    # todo failing. Odd because mean is passing; catastrophic cancellation in the online method?
    # other methods that use the variance estimator work just fine. Something about the gene issue
    # that is identified by other methods below?
    # (_cell_metrics, 'molecule_barcode_fraction_bases_above_30_variance',
    #  np.array(
    #      [np.nan, 0.0050, np.nan, np.nan, 0.0019, 0.0000, 0.0000, np.nan, 0.0015, np.nan, 0.0000,
    #       0.0000, np.nan, 0.0000, 0.0048, 0.0000, 0.0000, 0.0029, 0.0000, np.nan, 0.0000, 0.0044,
    #       0.0109, 0.0000, 0.0020, 0.0000, 0.0000, np.nan, 0.0000, 0.0100, np.nan, 0.0010, 0.0000,
    #       0.0052, 0.0000, 0.0000, 0.0000, 0.0011, 0.0000, 0.0162, 0.0000, 0.0016, 0.0000, np.nan,
    #       0.0178, 0.0020, np.nan, np.nan, 0.0000, 0.0163, np.nan, np.nan, 0.0011, np.nan, 0.0147,
    #       0.0018, 0.0007, 0.0306])),
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
    # todo right now the metrics count reads that have no 'gene' towards molecules, whereas
    # the calculations in the notebook exclude them. We should decide which method we prefer.
    # there may be further problems.
    # (_cell_metrics, 'reads_per_molecule',
    #  np.array(
    #      [1.0000, 2.0000, np.nan, 1.0000, 9.0000, 2.4000, 2.0000, 1.0000, 3.0000, 1.0000, 3.0000,
    #       3.0000, 1.0000, np.nan, 2.4167, 4.3333, 1.2222, 5.8750, 1.3333, 1.0000, 1.2000, 1.5000,
    #       4.6000, 2.0000, 2.5000, 1.2000, 2.1429, 1.0000, 2.6364, 4.0000, 1.0000, 2.1111, 1.7273,
    #       6.2500, 5.0000, 1.3333, 2.0000, 2.2500, np.nan, 2.0000, 4.3333, 3.9286, 2.2000, 1.0000,
    #       1.5000, 1.6667, np.nan, 1.0000, 1.6667, 1.8889, 1.0000, 1.0000, 2.2500, 1.0000, 9.7500,
    #       11.0000, 4.0000, 1.5000])),
    (_cell_metrics, 'reads_per_fragment',
     np.array(
         [1.0000, 1.0000, 1.0000, 1.0000, 1.1250, 1.3333, 2.0000, 1.0000, 1.2000, 1.0000, 1.2000,
          3.0000, 1.0000, 2.0000, 1.3182, 1.4444, 1.1000, 1.4688, 1.1429, 1.0000, 1.2000, 1.2857,
          1.5333, 2.0000, 1.2500, 1.0000, 1.1538, 1.0000, 1.3182, 1.0000, 1.0000, 1.4615, 1.3571,
          1.3158, 1.2500, 1.3333, 1.0000, 1.1250, 1.0000, 1.1765, 1.0833, 1.4103, 1.1000, 1.0000,
          1.2857, 1.2500, 1.0000, 1.0000, 1.2500, 1.3077, 1.0000, 1.0000, 1.2857, 1.0000, 1.3929,
          1.5714, 1.4737, 1.1053])),
    # (_cell_metrics, 'fragments_per_molecule',  # todo failure depends on above reads_per_molecule
    #  np.array(
    #      [1.0000, 2.0000, np.nan, 1.0000, 8.0000, 1.8000, 1.0000, 1.0000, 2.5000, 1.0000, 2.5000,
    #       1.0000, 1.0000, np.nan, 1.8333, 3.0000, 1.1111, 4.0000, 1.1667, 1.0000, 1.0000, 1.1667,
    #       3.0000, 1.0000, 2.0000, 1.2000, 1.8571, 1.0000, 2.0000, 4.0000, 1.0000, 1.4444, 1.2727,
    #       4.7500, 4.0000, 1.0000, 2.0000, 2.0000, np.nan, 1.7000, 4.0000, 2.7857, 2.0000, 1.0000,
    #       1.1667, 1.3333, np.nan, 1.0000, 1.3333, 1.4444, 1.0000, 1.0000, 1.7500, 1.0000, 7.0000,
    #       7.0000, 2.7143, 1.3571])),
    (_gene_metrics, 'molecule_barcode_fraction_bases_above_30_mean',
     np.array([1.0000, 1.0000, 0.8000, 0.9885, 0.9833, 0.9857, 0.7000, 0.9444])),
    (_gene_metrics, 'molecule_barcode_fraction_bases_above_30_variance',
     np.array([np.nan, np.nan, np.nan, 0.0011, 0.0051, 0.0014, np.nan, 0.0120])),
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
    """Test metrics that depend on other metrics

    This class tests a very large number of higher-order metrics that examine the functionality of
    the test suite across all measured instances of the metric class. E.g. for cell metrics (class),
    each test will verify the value for each cell (instance).

    Parameters
    ----------
    metrics : pd.DataFrame
        Output from subclass of sctools.metrics.MetricAggregator
    key : str
        The column of metrics to interrogate in the parametrized test
    expected_value : np.ndarray
        An array of expected values

    """
    # need to sort, metrics are not always in same order as results.
    observed = sorted(np.nan_to_num(metrics[key].values).round(4))
    expected_value = sorted(np.nan_to_num(expected_value))
    assert observed == expected_value


@pytest.mark.parametrize('metrics, key, expected_value', [
    # todo failing; suspect related to problem with how fragments are defined
    # (_cell_metrics, 'fragments_with_single_read_evidence', 345),
    # todo failing. Does not make sense that this would also be a fragment issue.
    # (_cell_metrics, 'molecules_with_single_read_evidence', 130),
    (_gene_metrics, 'fragments_with_single_read_evidence', 155),
    (_gene_metrics, 'molecules_with_single_read_evidence', 42),
])
def test_single_read_evidence(metrics, key, expected_value):
    """
    We want to determine how many molecules and fragments are covered by only one read, as reads
    covered by multiple reads have much lower probabilities of being the result of error processes.
    """
    observed = metrics[key].sum()
    assert observed == expected_value


def split_metrics_file(metrics_file):
    """
    produces two mergeable on-disk metric files from a single file that contain the first 3/4
    of the file in the first output and the last 3/4 of the file in the second output, such that
    1/2 of the metrics in the two files overlap
    """
    with fileinput.FileInput([metrics_file], mode='r', openhook=fileinput.hook_compressed) as f:
        data = [line for line in f]

    header, data = data[0], data[1:]

    low_split, high_split = round(len(data) * .25), round(len(data) * .75)
    file_1, file_2 = [_test_dir + 'metrics_for_merging_%d.csv' % i for i in (1, 2)]

    with open(file_1, 'wb') as f:
        f.write(header + b'\n')
        for line in data[:high_split]:
            f.write(line + b'\n')

    with open(file_2, 'wb') as f:
        f.write(header + b'\n')
        for line in data[low_split:]:
            f.write(line + b'\n')

    return file_1, file_2


@pytest.fixture
def mergeable_cell_metrics():
    return split_metrics_file(_cell_metric_output_file)


@pytest.fixture
def mergeable_gene_metrics():
    return split_metrics_file(_gene_metric_output_file)


def test_merge_cell_metrics_cli(mergeable_cell_metrics):
    """test the sctools merge cell metrics CLI invocation"""
    return_call = TenXV2.merge_cell_metrics(
        args=['-o', _test_dir + '/merged-cell-metrics.csv.gz'] + list(mergeable_cell_metrics))
    assert return_call == 0


def test_merge_gene_metrics_cli(mergeable_gene_metrics):
    """test the sctools merge gene metrics CLI invocation"""
    return_call = TenXV2.merge_gene_metrics(
        args=['-o', _test_dir + '/merged-gene-metrics.csv.gz'] + list(mergeable_gene_metrics))
    assert return_call == 0


def test_merge_cell_metrics_does_not_correct_duplicates(mergeable_cell_metrics):
    """
    test takes offset cell metrics outputs and merges them. Cell metrics does not check for
    duplication, so should return a 2x length file.
    """
    output_file = os.path.join(_test_dir, 'merged_metrics.csv.gz')
    m = MergeCellMetrics(mergeable_cell_metrics, output_file)
    m.execute()

    merged_data = pd.read_csv(output_file, index_col=0)

    input_sizes = []
    for f in mergeable_cell_metrics:
        input_sizes.append(pd.read_csv(f, index_col=0).shape)
    target_rows = sum(row for row, col in input_sizes)

    target_cols = input_sizes[0][1]  # cols will always be the same

    assert merged_data.shape == (target_rows, target_cols)


def test_merge_gene_metrics_averages_over_multiply_detected_genes(mergeable_gene_metrics):
    output_file = os.path.join(_test_dir, 'merged_metrics.csv.gz')
    m = MergeGeneMetrics(mergeable_gene_metrics, output_file)
    m.execute()

    merged_data = pd.read_csv(output_file, index_col=0)

    input_data = pd.read_csv(mergeable_gene_metrics[0], index_col=0)
    target_cols = input_data.shape[1]

    input_genes = input_data.index
    for f in mergeable_gene_metrics[1:]:
        input_genes = input_genes.union(pd.read_csv(f, index_col=0).index)
    target_rows = len(input_genes)

    assert merged_data.shape == (target_rows, target_cols), '%s' % repr(merged_data)


@pytest.mark.parametrize('bam, gatherer', [
    (_gene_sorted_bam, GatherGeneMetrics),
    (_cell_sorted_bam, GatherCellMetrics),
])
def test_gzip_compression(bam: str, gatherer: Callable):
    """
    gzip compression should produce a .gz file which is identical when uncompressed to the
    uncompressed version
    """

    gz_fout = _test_dir + 'test_bam.csv.gz'
    g: MetricGatherer = gatherer(bam, gz_fout, compress=True)
    g.extract_metrics()
    gz_metrics = pd.read_csv(gz_fout, index_col=0)

    fout = _test_dir + 'test_bam.csv'
    g: MetricGatherer = gatherer(bam, fout, compress=False)
    g.extract_metrics()
    metrics = pd.read_csv(fout, index_col=0)

    assert np.allclose(gz_metrics.fillna(0).values, metrics.fillna(0).values)
