import os
import pytest
import glob

import pysam

from .. import bam, platform


data_dir = os.path.split(__file__)[0] + '/data/'


# TEST SUBSETALIGNMENTS

@pytest.fixture(scope='module', params=[data_dir + 'test.sam', data_dir + 'test.bam'])
def sa_object(request):
    """fixture returns SubsetAlignments objects for testing"""
    return bam.SubsetAlignments(request.param)


@pytest.fixture(scope='module', params=[
    (data_dir + 'test.sam', 20, 20),
    (data_dir + 'test.bam', 20, 20)])
def indices(request):
    """fixture returns indices from a SubsetAlignments objects for testing"""
    sa = bam.SubsetAlignments(request.param[0])
    ind_specific, ind_nonspecific = sa.indices_by_chromosome(
        request.param[1], '19', include_other=request.param[2])
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
        n_specific, '19', include_other=n_nonspecific)
    assert len(ind_specific) == n_specific
    assert len(ind_nonspecific) == n_nonspecific


def test_indices_are_all_greater_than_zero(sa_object, n_specific, n_nonspecific):
    ind_specific, ind_nonspecific = sa_object.indices_by_chromosome(
        n_specific, '19', include_other=n_nonspecific)
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
    split_size = 0.02  # our test data is very small, 0.01mb = ~10kb, which should yield 5 files.
    with pytest.raises(RuntimeError):
        bam.split(bamfile, 'test_output', 'CB', approx_mb_per_split=split_size)


@pytest.fixture
def tagged_bam():
    args = [
        '--r1', data_dir + 'test_r1.fastq',
        '--i1', data_dir + 'test_i7.fastq',
        '--u2', data_dir + 'test_r2.bam',
        '--output-bamfile', 'test_tagged_bam.bam',
        '--whitelist', data_dir + '1k-august-2016.txt']
    platform.TenXV2.attach_barcodes(args)
    return 'test_tagged_bam.bam'


def test_split_on_tagged_bam(tagged_bam):
    split_size = 0.005  # our test data is very small, this value should yield 3 files
    outputs = bam.split(tagged_bam, 'test_output', 'CB', 'CR', approx_mb_per_split=split_size)
    assert len(outputs) == 3

    # cleanup
    os.remove(tagged_bam)  # clean up
    for f in glob.glob('test_output_*'):
        os.remove(f)


def test_split_with_large_chunk_size_generates_one_file(tagged_bam):
    split_size = 1024  # our test data is very small, this value should yield 1 file
    outputs = bam.split(tagged_bam, 'test_output', 'CB', 'CR', approx_mb_per_split=split_size)
    assert len(outputs) == 1

    # the file should be full size
    with pysam.AlignmentFile(outputs[0], 'rb', check_sq=False) as f:
        assert len([x for x in f]) == 100

    # cleanup
    os.remove(tagged_bam)  # clean up
    for f in glob.glob('test_output_*'):
        os.remove(f)


def test_split_with_raise_missing_true_raises_warning_without_cr_barcode_passed(tagged_bam):
    split_size = 1024  # our test data is very small, this value should yield 1 file
    with pytest.raises(RuntimeError):
        outputs = bam.split(tagged_bam, 'test_output', 'CB', approx_mb_per_split=split_size,
                            raise_missing=True)

    # cleanup
    os.remove(tagged_bam)  # clean up
    for f in glob.glob('test_output_*'):
        os.remove(f)


def test_split_succeeds_with_raise_missing_false_and_no_cr_barcode_passed(tagged_bam):
    split_size = 1024  # our test data is very small, this value should yield 1 file
    outputs = bam.split(tagged_bam, 'test_output', 'CB', approx_mb_per_split=split_size,
                        raise_missing=False)

    assert len(outputs) == 1

    # the file should be full size
    with pysam.AlignmentFile(outputs[0], 'rb', check_sq=False) as f:
        assert len([x for x in f]) == 1  # only one of our barcodes is whitelisted or within 1 base

    # cleanup
    os.remove(tagged_bam)  # clean up
    for f in glob.glob('test_output_*'):
        os.remove(f)
