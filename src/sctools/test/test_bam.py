import pytest
import os
from .. import bam


data_dir = os.path.split(__file__)[0] + '/data'


# TEST SUBSETALIGNMENTS

@pytest.fixture(scope='module', params=[data_dir + '/test.sam', data_dir + '/test.bam'])
def sa_object(request):
    """fixture returns SubsetAlignments objects for testing"""
    return bam.SubsetAlignments(request.param)


@pytest.fixture(scope='module', params=[
    (data_dir + '/test.sam', 20, 20),
    (data_dir + '/test.bam', 20, 20)])
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
