from copy import deepcopy
import glob
import os
import shutil

import pysam
import pytest

from .. import bam, platform, consts

data_dir = os.path.split(__file__)[0] + '/data/'


TAG_KEYS = ['FOO', 'BAR', 'BAZ']

SORTED_VALUES = [
    (['A', 'A', 'A'], 'A'),
    (['A', 'A', 'A'], 'B'),
    (['A', 'A', 'B'], 'A'),
    (['A', 'A', 'B'], 'B'),
    (['A', 'B', 'A'], 'A'),
    (['A', 'B', 'A'], 'B'),
    (['A', 'B', 'B'], 'A'),
    (['A', 'B', 'B'], 'B'),
    (['B', 'A', 'A'], 'A'),
    (['B', 'A', 'A'], 'B'),
    (['B', 'A', 'B'], 'A'),
    (['B', 'A', 'B'], 'B'),
    (['B', 'B', 'A'], 'A'),
    (['B', 'B', 'A'], 'B'),
    (['B', 'B', 'B'], 'A'),
    (['B', 'B', 'B'], 'B'),
]

UNSORTED_VALUES = [
    (['B', 'B', 'A'], 'A'),
    (['B', 'B', 'B'], 'A'),
    (['B', 'A', 'A'], 'A'),
    (['B', 'B', 'A'], 'B'),
    (['B', 'B', 'B'], 'B'),
    (['A', 'A', 'B'], 'A'),
    (['A', 'B', 'A'], 'B'),
    (['A', 'B', 'B'], 'A'),
    (['A', 'A', 'A'], 'B'),
    (['B', 'A', 'B'], 'A'),
    (['A', 'B', 'A'], 'A'),
    (['B', 'A', 'A'], 'B'),
    (['A', 'A', 'A'], 'A'),
    (['A', 'A', 'B'], 'B'),
    (['A', 'B', 'B'], 'B'),
    (['B', 'A', 'B'], 'B'),
]


# TEST SUBSETALIGNMENTS


@pytest.fixture(scope='module', params=[data_dir + 'test.sam', data_dir + 'test.bam'])
def sa_object(request):
    """fixture returns SubsetAlignments objects for testing"""
    return bam.SubsetAlignments(request.param)


@pytest.fixture(
    scope='module',
    params=[(data_dir + 'test.sam', 20, 20), (data_dir + 'test.bam', 20, 20)],
)
def indices(request):
    """fixture returns indices from a SubsetAlignments objects for testing"""
    sa = bam.SubsetAlignments(request.param[0])
    ind_specific, ind_nonspecific = sa.indices_by_chromosome(
        request.param[1], '19', include_other=request.param[2]
    )
    return ind_specific, ind_nonspecific


@pytest.fixture(scope='session')
def n_nonspecific():
    """the number of non-specific records to extract"""
    return 20


@pytest.fixture(scope='session')
def n_specific():
    """the number of specific records to extract"""
    return 20


def test_incorrect_extension_does_not_raise_when_open_mode_is_specified():
    sa = bam.SubsetAlignments(data_dir + '/is_actually_a_samfile.wrong_extensions', 'r')
    assert isinstance(sa, bam.SubsetAlignments)


def test_incorrect_extension_without_open_mode_raises_value_error():
    with pytest.raises(ValueError):
        _ = bam.SubsetAlignments(data_dir + '/is_actually_a_samfile.wrong_extensions')


def test_str_and_int_chromosomes_both_function(sa_object):
    _ = sa_object.indices_by_chromosome(10, '19', include_other=10)
    _ = sa_object.indices_by_chromosome(10, 19, include_other=10)


def test_correct_number_of_indices_are_extracted(sa_object, n_specific, n_nonspecific):
    ind_specific, ind_nonspecific = sa_object.indices_by_chromosome(
        n_specific, '19', include_other=n_nonspecific
    )
    assert len(ind_specific) == n_specific
    assert len(ind_nonspecific) == n_nonspecific


def test_indices_are_all_greater_than_zero(sa_object, n_specific, n_nonspecific):
    ind_specific, ind_nonspecific = sa_object.indices_by_chromosome(
        n_specific, '19', include_other=n_nonspecific
    )
    assert min(ind_nonspecific + ind_nonspecific) > 0


def test_chromosome_19_comes_before_21(indices):
    """chromosome 19 comes before 21 in the test file, this should be replicated in the output"""
    assert max(indices[0]) < min(indices[1])


# TAGGER TESTED IN INTEGRATION TESTS ONLY (see test_entrypoints.py)

# TEST SPLIT


@pytest.fixture(scope='module', params=[data_dir + 'test.sam', data_dir + 'test.bam'])
def bamfile(request):
    return request.param


def test_split_bam_raises_value_error_when_passed_bam_without_barcodes(bamfile):
    split_size = (
        0.02
    )  # our test data is very small, 0.01mb = ~10kb, which should yield 5 files.
    with pytest.raises(RuntimeError):
        bam.split(
            [bamfile],
            'test_output',
            [consts.CELL_BARCODE_TAG_KEY],
            approx_mb_per_split=split_size,
        )


@pytest.fixture
def tagged_bam():
    args = [
        '--r1',
        data_dir + 'test_r1.fastq',
        '--i1',
        data_dir + 'test_i7.fastq',
        '--u2',
        data_dir + 'test_r2.bam',
        '--output-bamfile',
        'test_tagged_bam.bam',
        '--whitelist',
        data_dir + '1k-august-2016.txt',
    ]
    platform.TenXV2.attach_barcodes(args)
    return 'test_tagged_bam.bam'


def test_split_on_tagged_bam(tagged_bam):
    split_size = 0.005  # our test data is very small, this value should yield 3 files
    outputs = bam.split(
        [tagged_bam],
        'test_output',
        [consts.CELL_BARCODE_TAG_KEY, consts.RAW_CELL_BARCODE_TAG_KEY],
        approx_mb_per_split=split_size,
    )
    assert len(outputs) == 3

    # cleanup
    os.remove(tagged_bam)  # clean up
    for f in glob.glob('test_output_*'):
        os.remove(f)


def test_split_with_large_chunk_size_generates_one_file(tagged_bam):
    split_size = 1024  # our test data is very small, this value should yield 1 file
    outputs = bam.split(
        [tagged_bam],
        'test_output',
        [consts.CELL_BARCODE_TAG_KEY, consts.RAW_CELL_BARCODE_TAG_KEY],
        approx_mb_per_split=split_size,
    )
    assert len(outputs) == 1

    # the file should be full size
    with pysam.AlignmentFile(outputs[0], 'rb', check_sq=False) as f:
        assert len([x for x in f]) == 100

    # cleanup
    os.remove(tagged_bam)  # clean up
    for f in glob.glob('test_output_*'):
        os.remove(f)


def test_split_with_raise_missing_true_raises_warning_without_cr_barcode_passed(
    tagged_bam
):
    split_size = 1024  # our test data is very small, this value should yield 1 file
    with pytest.raises(RuntimeError):
        bam.split(
            [tagged_bam],
            'test_output',
            [consts.CELL_BARCODE_TAG_KEY],
            approx_mb_per_split=split_size,
            raise_missing=True,
        )

    # cleanup
    os.remove(tagged_bam)  # clean up
    for f in glob.glob('test_output_*'):
        os.remove(f)


def test_split_succeeds_with_raise_missing_false_and_no_cr_barcode_passed(tagged_bam):
    split_size = 1024  # our test data is very small, this value should yield 1 file
    outputs = bam.split(
        [tagged_bam],
        'test_output',
        [consts.CELL_BARCODE_TAG_KEY],
        approx_mb_per_split=split_size,
        raise_missing=False,
    )

    assert len(outputs) == 1

    # the file should be full size
    with pysam.AlignmentFile(outputs[0], 'rb', check_sq=False) as f:
        assert (
            len([x for x in f]) == 1
        )  # only one of our barcodes is whitelisted or within 1 base

    # cleanup
    os.remove(tagged_bam)  # clean up
    for f in glob.glob('test_output_*'):
        os.remove(f)


def test_get_barcodes_from_bam(tagged_bam):
    outputs = bam.get_barcodes_from_bam(
        tagged_bam,
        [consts.CELL_BARCODE_TAG_KEY, consts.RAW_CELL_BARCODE_TAG_KEY],
        raise_missing=True,
    )
    assert len(outputs) == 99


def test_get_barcodes_from_bam_with_raise_missing_true_raises_warning_without_cr_barcode_passed(
    tagged_bam
):
    with pytest.raises(RuntimeError):
        bam.get_barcodes_from_bam(
            tagged_bam, [consts.CELL_BARCODE_TAG_KEY], raise_missing=True
        )


def test_write_barcodes_to_bins(tagged_bam):
    barcodes = bam.get_barcodes_from_bam(
        tagged_bam,
        [consts.CELL_BARCODE_TAG_KEY, consts.RAW_CELL_BARCODE_TAG_KEY],
        raise_missing=True,
    )

    test_barcodes_to_bins = {}
    for barcode in barcodes:
        test_barcodes_to_bins[barcode] = 0

    filenames = bam.write_barcodes_to_bins(
        tagged_bam,
        [consts.CELL_BARCODE_TAG_KEY, consts.RAW_CELL_BARCODE_TAG_KEY],
        test_barcodes_to_bins,
        raise_missing=False,
    )

    assert len(filenames) == 1

    # cleanup
    for f in filenames:
        shutil.rmtree(os.path.dirname(f))


def test_get_barcode_for_alignment(tagged_bam):
    with pysam.AlignmentFile(tagged_bam, 'rb', check_sq=False) as input_alignments:
        for alignment in input_alignments:
            barcode = bam.get_barcode_for_alignment(
                alignment,
                [consts.CELL_BARCODE_TAG_KEY, consts.RAW_CELL_BARCODE_TAG_KEY],
                raise_missing=False,
            )
            assert barcode == "NTAAGAGTCTGCAAGT"
            break


def test_get_barcode_for_alignment_raises_error_for_missing_tag(tagged_bam):
    with pysam.AlignmentFile(tagged_bam, 'rb', check_sq=False) as input_alignments:
        for alignment in input_alignments:
            with pytest.raises(RuntimeError):
                bam.get_barcode_for_alignment(alignment, TAG_KEYS, raise_missing=True)


# TEST SORTING


def test_tag_sortable_records_compare_correctly():
    records = make_records_from_values(TAG_KEYS, SORTED_VALUES)
    num_records = len(SORTED_VALUES)
    for i in range(num_records):
        for j in range(num_records):
            if i < j:
                assert records[i] < records[j]
            elif i == j:
                assert records[i] == records[j]
            else:
                assert records[i] > records[j]


def test_tag_sortable_records_raises_error_on_different_tag_lists():
    r1 = bam.TagSortableRecord(['FOO', 'BAR'], ['A', 'A'], 'A')
    r2 = bam.TagSortableRecord(['BAR', 'BAZ'], ['A', 'A'], 'A')
    with pytest.raises(ValueError):
        r1 == r2


def test_tag_sortable_records_str():
    record = bam.TagSortableRecord(TAG_KEYS, SORTED_VALUES[0][0], SORTED_VALUES[0][1])
    s = record.__str__()
    assert 'TagSortableRecord' in s
    assert "['FOO', 'BAR', 'BAZ']" in s


def test_verify_sort_on_unsorted_records_raises_error():
    records = make_records_from_values(TAG_KEYS, UNSORTED_VALUES)
    with pytest.raises(bam.SortError):
        bam.verify_sort(records, TAG_KEYS)


def test_verify_sort_raises_no_error_on_sorted_records():
    records = make_records_from_values(TAG_KEYS, SORTED_VALUES)
    bam.verify_sort(records, TAG_KEYS)


def test_sort_by_tags_and_queryname_sorts_correctly_from_file():
    tag_keys = ['UB', 'CB', 'GE']
    with pysam.AlignmentFile(data_dir + 'unsorted.bam', 'rb') as f:
        records = f.fetch(until_eof=True)
        sorted_records = bam.sort_by_tags_and_queryname(records, tag_keys)
    tag_sortable_records = (
        bam.TagSortableRecord.from_aligned_segment(r, tag_keys) for r in sorted_records
    )
    bam.verify_sort(tag_sortable_records, tag_keys)


def test_sort_by_tags_and_queryname_sorts_correctly_from_file_no_tag_keys():
    tag_keys = []
    with pysam.AlignmentFile(data_dir + 'unsorted.bam', 'rb') as f:
        records = f.fetch(until_eof=True)
        sorted_records = bam.sort_by_tags_and_queryname(records, tag_keys)
    tag_sortable_records = (
        bam.TagSortableRecord.from_aligned_segment(r, tag_keys) for r in sorted_records
    )
    bam.verify_sort(tag_sortable_records, tag_keys)


def test_tag_sortable_records_sort_correctly():
    tag_keys = TAG_KEYS
    records = make_records_from_values(tag_keys, deepcopy(UNSORTED_VALUES))
    sorted_records = sorted(records)
    bam.verify_sort(sorted_records, tag_keys)


def test_tag_sortable_records_sort_correctly_when_already_sorted():
    # This is to a bit paranoid, but just make sure sorted stays correct if already sorted
    tag_keys = TAG_KEYS
    records = make_records_from_values(tag_keys, deepcopy(SORTED_VALUES))
    sorted_records = sorted(records)
    bam.verify_sort(sorted_records, tag_keys)


def test_sort_by_tags_and_queryname_sorts_correctly_no_tag_keys():
    tag_keys = []
    records = make_records_from_values(tag_keys, deepcopy(UNSORTED_VALUES))
    sorted_records = sorted(records)
    bam.verify_sort(sorted_records, tag_keys)


def test_tag_sortable_record_missing_tag_value_is_empty_string():
    tags = ['_NOT_REAL_TAG_']
    with pysam.AlignmentFile(data_dir + 'unsorted.bam', 'rb') as f:
        records = f.fetch(until_eof=True)
        first_record = next(iter(records))
        sortable_record = bam.TagSortableRecord.from_aligned_segment(first_record, tags)
        assert sortable_record.tag_values[0] == ''


def test_tag_sortable_record_lt_is_false_for_equal_records():
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    assert not r1 < r2


def test_tag_sortable_record_lt_is_true_for_smaller_query_name():
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='B'
    )
    assert r1 < r2


def test_tag_sortable_record_lt_is_true_for_smaller_tag():
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'B'], query_name='A'
    )
    assert r1 < r2
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'B', 'A'], query_name='A'
    )
    assert r1 < r2
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['B', 'A', 'A'], query_name='A'
    )
    assert r1 < r2


def test_tag_sortable_record_lt_is_true_for_smaller_tag_regardless_of_query_name():
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='B'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'B'], query_name='A'
    )
    assert r1 < r2


def test_tag_sortable_record_lt_empty_query_name_is_smaller():
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name=''
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    assert r1 < r2


def test_tag_sortable_record_lt_empty_tag_is_smaller():
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', ''], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    assert r1 < r2
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', '', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    assert r1 < r2
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['', 'A', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    assert r1 < r2


def test_tag_sortable_record_eq_is_true_for_identical_records():
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    assert r1 == r2


def test_tag_sortable_record_eq_is_false_when_any_difference_exists():
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='B'
    )
    assert not r1 == r2
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'B'], query_name='A'
    )
    assert not r1 == r2
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'B', 'A'], query_name='A'
    )
    assert not r1 == r2
    r1 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['A', 'A', 'A'], query_name='A'
    )
    r2 = bam.TagSortableRecord(
        tag_keys=TAG_KEYS, tag_values=['B', 'A', 'A'], query_name='A'
    )
    assert not r1 == r2


def make_records_from_values(tag_keys, tags_and_query_name):
    records = []
    for i in range(len(tags_and_query_name)):
        r = bam.TagSortableRecord(
            tag_keys=tag_keys,
            tag_values=tags_and_query_name[i][0],
            query_name=tags_and_query_name[i][1],
        )
        records.append(r)
    return records
