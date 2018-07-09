import glob
import os
import tempfile

import numpy as np
import pysam
import scipy.sparse as sp

from sctools import platform, count, consts

data_dir = os.path.split(__file__)[0] + '/data/'


def test_Attach10XBarcodes_entrypoint():
    args = [
        '--r1', data_dir + 'test_r1.fastq',
        '--i1', data_dir + 'test_i7.fastq',
        '--u2', data_dir + 'test.bam',
        '--output-bamfile', 'test_tagged_bam.bam']

    rc = platform.TenXV2.attach_barcodes(args)
    assert rc == 0
    with pysam.AlignmentFile('test_tagged_bam.bam', 'rb', check_sq=False) as f:
        for alignment in f:
            # each alignment should now have a tag, and that tag should be a string
            assert isinstance(alignment.get_tag(consts.QUALITY_CELL_BARCODE_TAG_KEY), str)
            assert isinstance(alignment.get_tag(consts.RAW_CELL_BARCODE_TAG_KEY), str)
            assert isinstance(alignment.get_tag(consts.QUALITY_MOLECULE_BARCODE_TAG_KEY), str)
            assert isinstance(alignment.get_tag(consts.RAW_MOLECULE_BARCODE_TAG_KEY), str)
            assert isinstance(alignment.get_tag(consts.RAW_SAMPLE_BARCODE_TAG_KEY), str)
            assert isinstance(alignment.get_tag(consts.QUALITY_SAMPLE_BARCODE_TAG_KEY), str)
    os.remove('test_tagged_bam.bam')  # clean up


def test_Attach10XBarcodes_entrypoint_with_whitelist():
    args = [
        '--r1', data_dir + 'test_r1.fastq',
        '--i1', data_dir + 'test_i7.fastq',
        '--u2', data_dir + 'test.bam',
        '--output-bamfile', 'test_tagged_bam.bam',
        '--whitelist', data_dir + '1k-august-2016.txt']

    return_call = platform.TenXV2.attach_barcodes(args)
    assert return_call == 0
    success = False
    with pysam.AlignmentFile('test_tagged_bam.bam', 'rb', check_sq=False) as f:
        for alignment in f:
            if alignment.has_tag(consts.CELL_BARCODE_TAG_KEY):
                success = True
            # each alignment should now have a tag, and that tag should be a string
            assert isinstance(alignment.get_tag(consts.RAW_CELL_BARCODE_TAG_KEY), str)
            assert isinstance(alignment.get_tag(consts.QUALITY_CELL_BARCODE_TAG_KEY), str)
            assert isinstance(alignment.get_tag(consts.RAW_MOLECULE_BARCODE_TAG_KEY), str)
            assert isinstance(alignment.get_tag(consts.QUALITY_MOLECULE_BARCODE_TAG_KEY), str)
            assert isinstance(alignment.get_tag(consts.RAW_SAMPLE_BARCODE_TAG_KEY), str)
            assert isinstance(alignment.get_tag(consts.QUALITY_SAMPLE_BARCODE_TAG_KEY), str)
    assert success
    os.remove('test_tagged_bam.bam')  # clean up


def test_split_bam():
    tag_args = [
        '--r1', data_dir + 'test_r1.fastq',
        '--i1', data_dir + 'test_i7.fastq',
        '--u2', data_dir + 'test.bam',
        '--output-bamfile', 'test_tagged_bam.bam',
        '--whitelist', data_dir + '1k-august-2016.txt']

    platform.TenXV2.attach_barcodes(tag_args)

    split_args = [
        '--bamfile', 'test_tagged_bam.bam',
        '--output-prefix', 'test_tagged',
        '--subfile-size', '0.005',
        '--tags', consts.CELL_BARCODE_TAG_KEY, consts.RAW_CELL_BARCODE_TAG_KEY]

    return_call = platform.GenericPlatform.split_bam(split_args)
    assert return_call == 0

    for f in glob.glob('test_tagged*'):
        os.remove(f)


def test_count_merge():
    tmp = tempfile.mkdtemp()

    data, ind, col = [np.arange(10)] * 3
    matrix = sp.coo_matrix((data, (ind, col)), shape=(10, 10), dtype=np.float32).tocsr()
    # be lazy and reuse the inds as the col and row index
    counts = count.CountMatrix(matrix, ind, col)
    counts.save(tmp + '/test_input_1')
    counts.save(tmp + '/test_input_2')

    merge_args = [
        '-o', tmp + '/test_merged_counts',
        '-i', tmp + '/test_input_2', tmp + '/test_input_1'
    ]
    return_call = platform.GenericPlatform.merge_count_matrices(merge_args)
    assert return_call == 0
