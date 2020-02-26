import glob
import os
import tempfile

import numpy as np
import pysam
import pytest
import scipy.sparse as sp

from sctools import bam, platform, count, consts

data_dir = os.path.split(__file__)[0] + "/data/"


def test_Attach10XBarcodes_entrypoint():
    args = [
        "--r1",
        data_dir + "test_r1.fastq",
        "--i1",
        data_dir + "test_i7.fastq",
        "--u2",
        data_dir + "test.bam",
        "--output-bamfile",
        "test_tagged_bam.bam",
    ]

    rc = platform.TenXV2.attach_barcodes(args)
    assert rc == 0
    with pysam.AlignmentFile("test_tagged_bam.bam", "rb", check_sq=False) as f:
        for alignment in f:
            # each alignment should now have a tag, and that tag should be a string
            assert isinstance(
                alignment.get_tag(consts.QUALITY_CELL_BARCODE_TAG_KEY), str
            )
            assert isinstance(alignment.get_tag(consts.RAW_CELL_BARCODE_TAG_KEY), str)
            assert isinstance(
                alignment.get_tag(consts.QUALITY_MOLECULE_BARCODE_TAG_KEY), str
            )
            assert isinstance(
                alignment.get_tag(consts.RAW_MOLECULE_BARCODE_TAG_KEY), str
            )
            assert isinstance(alignment.get_tag(consts.RAW_SAMPLE_BARCODE_TAG_KEY), str)
            assert isinstance(
                alignment.get_tag(consts.QUALITY_SAMPLE_BARCODE_TAG_KEY), str
            )
    os.remove("test_tagged_bam.bam")  # clean up


def test_Attach10XBarcodes_entrypoint_with_whitelist():
    args = [
        "--r1",
        data_dir + "test_r1.fastq",
        "--i1",
        data_dir + "test_i7.fastq",
        "--u2",
        data_dir + "test.bam",
        "--output-bamfile",
        "test_tagged_bam.bam",
        "--whitelist",
        data_dir + "1k-august-2016.txt",
    ]

    return_call = platform.TenXV2.attach_barcodes(args)
    assert return_call == 0
    success = False
    with pysam.AlignmentFile("test_tagged_bam.bam", "rb", check_sq=False) as f:
        for alignment in f:
            if alignment.has_tag(consts.CELL_BARCODE_TAG_KEY):
                success = True
            # each alignment should now have a tag, and that tag should be a string
            assert isinstance(alignment.get_tag(consts.RAW_CELL_BARCODE_TAG_KEY), str)
            assert isinstance(
                alignment.get_tag(consts.QUALITY_CELL_BARCODE_TAG_KEY), str
            )
            assert isinstance(
                alignment.get_tag(consts.RAW_MOLECULE_BARCODE_TAG_KEY), str
            )
            assert isinstance(
                alignment.get_tag(consts.QUALITY_MOLECULE_BARCODE_TAG_KEY), str
            )
            assert isinstance(alignment.get_tag(consts.RAW_SAMPLE_BARCODE_TAG_KEY), str)
            assert isinstance(
                alignment.get_tag(consts.QUALITY_SAMPLE_BARCODE_TAG_KEY), str
            )
    assert success
    os.remove("test_tagged_bam.bam")  # clean up


def test_AttachBarcodes_entrypoint_with_whitelist():
    # test of the BarcodePlatform.attach_barcodes entry point with
    # sample, cell, and molecule barcodes all specified
    args = [
        "--r1",
        data_dir + "test_r1.fastq",
        "--i1",
        data_dir + "test_i7.fastq",
        "--u2",
        data_dir + "test.bam",
        "--output-bamfile",
        "test_tagged_bam.bam",
        "--whitelist",
        data_dir + "1k-august-2016.txt",
        "--sample-barcode-start-position",
        "0",
        "--sample-barcode-length",
        "8",
        "--cell-barcode-start-position",
        "0",
        "--cell-barcode-length",
        "16",
        "--molecule-barcode-start-position",
        "16",
        "--molecule-barcode-length",
        "7",  # changed 10>7 intentionally for test
    ]

    return_call = platform.BarcodePlatform.attach_barcodes(args)
    assert return_call == 0
    success = False
    with pysam.AlignmentFile("test_tagged_bam.bam", "rb", check_sq=False) as f:
        for alignment in f:
            if alignment.has_tag(consts.CELL_BARCODE_TAG_KEY):
                success = True
            # each alignment should now have a tag, and that tag should be a string
            assert isinstance(alignment.get_tag(consts.RAW_CELL_BARCODE_TAG_KEY), str)
            assert isinstance(
                alignment.get_tag(consts.QUALITY_CELL_BARCODE_TAG_KEY), str
            )
            assert isinstance(
                alignment.get_tag(consts.RAW_MOLECULE_BARCODE_TAG_KEY), str
            )
            assert len(alignment.get_tag(consts.RAW_MOLECULE_BARCODE_TAG_KEY)) == 7
            assert isinstance(
                alignment.get_tag(consts.QUALITY_MOLECULE_BARCODE_TAG_KEY), str
            )
            assert isinstance(alignment.get_tag(consts.RAW_SAMPLE_BARCODE_TAG_KEY), str)
            assert isinstance(
                alignment.get_tag(consts.QUALITY_SAMPLE_BARCODE_TAG_KEY), str
            )
    assert success
    os.remove("test_tagged_bam.bam")  # clean up


def test_split_bam():
    tag_args = [
        "--r1",
        data_dir + "test_r1.fastq",
        "--i1",
        data_dir + "test_i7.fastq",
        "--u2",
        data_dir + "test.bam",
        "--output-bamfile",
        "test_tagged_bam.bam",
        "--whitelist",
        data_dir + "1k-august-2016.txt",
    ]

    platform.TenXV2.attach_barcodes(tag_args)

    split_args = [
        "--bamfile",
        "test_tagged_bam.bam",
        "--output-prefix",
        "test_tagged",
        "--subfile-size",
        "0.005",
        "--tags",
        consts.CELL_BARCODE_TAG_KEY,
        consts.RAW_CELL_BARCODE_TAG_KEY,
    ]

    return_call = platform.GenericPlatform.split_bam(split_args)
    assert return_call == 0

    for f in glob.glob("test_tagged*"):
        os.remove(f)


def test_tag_sort_bam():
    args = [
        "-i",
        data_dir + "unsorted.bam",
        "-o",
        "test_sorted.bam",
        "-t",
        consts.CELL_BARCODE_TAG_KEY,
        consts.GENE_NAME_TAG_KEY,
        consts.MOLECULE_BARCODE_TAG_KEY,
    ]

    return_call = platform.GenericPlatform.tag_sort_bam(args)
    assert return_call == 0

    tag_keys = [
        consts.CELL_BARCODE_TAG_KEY,
        consts.GENE_NAME_TAG_KEY,
        consts.MOLECULE_BARCODE_TAG_KEY,
    ]
    with pysam.AlignmentFile("test_sorted.bam", "rb") as f:
        segments = f.fetch(until_eof=True)
        tag_sortable_records = (
            bam.TagSortableRecord.from_aligned_segment(s, tag_keys) for s in segments
        )
        bam.verify_sort(tag_sortable_records, tag_keys)

    for f in glob.glob("test_sorted*"):
        os.remove(f)


def test_tag_sort_bam_dash_t_specified_multiple_times():
    args = [
        "-i",
        data_dir + "unsorted.bam",
        "-o",
        "test_sorted.bam",
        "-t",
        consts.CELL_BARCODE_TAG_KEY,
        "-t",
        consts.GENE_NAME_TAG_KEY,
        "-t",
        consts.MOLECULE_BARCODE_TAG_KEY,
    ]

    return_call = platform.GenericPlatform.tag_sort_bam(args)
    assert return_call == 0

    tag_keys = [
        consts.CELL_BARCODE_TAG_KEY,
        consts.GENE_NAME_TAG_KEY,
        consts.MOLECULE_BARCODE_TAG_KEY,
    ]
    with pysam.AlignmentFile("test_sorted.bam", "rb") as f:
        segments = f.fetch(until_eof=True)
        tag_sortable_record_generator = (
            bam.TagSortableRecord.from_aligned_segment(s, tag_keys) for s in segments
        )
        bam.verify_sort(tag_sortable_record_generator, tag_keys)

    for f in glob.glob("test_sorted*"):
        os.remove(f)


def test_tag_sort_bam_no_tags():
    args = ["-i", data_dir + "unsorted.bam", "-o", "test_sorted.bam"]

    return_call = platform.GenericPlatform.tag_sort_bam(args)
    assert return_call == 0

    tag_keys = []
    with pysam.AlignmentFile("test_sorted.bam", "rb") as f:
        segments = f.fetch(until_eof=True)
        tag_sortable_records = (
            bam.TagSortableRecord.from_aligned_segment(s, tag_keys) for s in segments
        )
        bam.verify_sort(tag_sortable_records, tag_keys)

    for f in glob.glob("test_sorted*"):
        os.remove(f)


def test_verify_bam_sort():
    args = [
        "-i",
        data_dir + "cell-gene-umi-queryname-sorted.bam",
        "-t",
        consts.CELL_BARCODE_TAG_KEY,
        consts.GENE_NAME_TAG_KEY,
        consts.MOLECULE_BARCODE_TAG_KEY,
    ]

    return_call = platform.GenericPlatform.verify_bam_sort(args)
    assert return_call == 0


def test_verify_bam_sort_raises_error_on_unsorted():
    args = [
        "-i",
        data_dir + "unsorted.bam",
        "-t",
        consts.CELL_BARCODE_TAG_KEY,
        consts.GENE_NAME_TAG_KEY,
        consts.MOLECULE_BARCODE_TAG_KEY,
    ]

    with pytest.raises(bam.SortError):
        platform.GenericPlatform.verify_bam_sort(args)


def test_count_merge():
    tmp = tempfile.mkdtemp()

    data, ind, col = [np.arange(10)] * 3
    matrix = sp.coo_matrix((data, (ind, col)), shape=(10, 10), dtype=np.float32).tocsr()
    # be lazy and reuse the inds as the col and row index
    counts = count.CountMatrix(matrix, ind, col)
    counts.save(tmp + "/test_input_1")
    counts.save(tmp + "/test_input_2")

    merge_args = [
        "-o",
        tmp + "/test_merged_counts",
        "-i",
        tmp + "/test_input_2",
        tmp + "/test_input_1",
    ]
    return_call = platform.GenericPlatform.merge_count_matrices(merge_args)
    assert return_call == 0
