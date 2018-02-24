import os
import pytest
import tempfile
from sctools import platform

_data_dir = os.path.split(__file__)[0] + '/data'
_gene_sorted_bam = _data_dir + '/gene-sorted.bam'
_cell_sorted_bam = _data_dir + '/cell-sorted.bam'

_test_dir = tempfile.mkdtemp()


def test_calculate_cell_metrics():
    platform.TenXV2.calculate_cell_metrics(
        args=['-i', _cell_sorted_bam, '-o', _test_dir + '/gene_metrics.csv'])


def test_calculate_gene_metrics():
    platform.TenXV2.calculate_gene_metrics(
        args=['-i', _gene_sorted_bam, '-o', _test_dir + '/gene_metrics.csv'])


def test_merge_cell_metrics():
    platform.TenXV2.merge_cell_metrics(
        args=['-o', _test_dir + '/merged-cell-metrics.csv',
              'cell_metrics_1.csv', 'cell_metrics_2.csv']
    )


def test_merge_gene_metrics():
    platform.TenXV2.merge_gene_metrics(
        args=['-o', _test_dir + '/merged-gene-metrics.csv',
              'gene_metrics_1.csv', 'gene_metrics_2.csv']
    )
