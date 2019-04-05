import os
from .. import gtf
from itertools import chain
import pytest

_data_dir = os.path.split(__file__)[0] + '/data'
_files = ['%s/%s' % (_data_dir, f) for f in ('test.gtf', 'test.gtf.gz', 'test.gtf.bz2')]


@pytest.fixture(scope='module', params=_files)
def files(request):
    """returns a filename"""
    return request.param


def test_opens_file_reads_first_line(files):
    rd = gtf.Reader(files, 'r', header_comment_char='#')
    record = next(iter(rd))
    assert isinstance(record, gtf.GTFRecord)


def test_opens_file_populates_fields_properly(files):
    rd = gtf.Reader(files, 'r', header_comment_char='#')
    record = next(iter(rd))
    assert record.seqname == 'chr19'
    assert record.chromosome == 'chr19'
    assert record.source == 'HAVANA'
    assert record.feature == 'gene'
    assert record.start == 60951
    assert record.end == 71626
    assert record.score == '.'
    assert record.strand == '-'
    assert record.frame == '.'

    expected_features = {
        'gene_id': 'ENSG00000282458.1',
        'gene_type': 'transcribed_processed_pseudogene',
        'gene_status': 'KNOWN',
        'gene_name': 'WASH5P',
        'level': '2',
        'havana_gene': 'OTTHUMG00000180466.8',
    }
    assert record._attributes == expected_features

    assert all(
        i in str(record)
        for i in chain(expected_features.keys(), expected_features.values())
    )


def test_set_attribute_verify_included_in_output_string(files):
    rd = gtf.Reader(files, 'r', header_comment_char='#')
    record = next(iter(rd))
    record.set_attribute('test_attr', 'foo')
    assert record.get_attribute('test_attr') == 'foo'

    # verify in output string
    assert 'foo' in str(record)


def test_opens_file_parses_size(files):
    rd = gtf.Reader(files, 'r', header_comment_char='#')
    record = next(iter(rd))
    assert 71626 - 60951 == record.size

    # mangle record, make sure error is raised
    record._fields[3:5] = [record.end, record.start]
    with pytest.raises(ValueError):
        getattr(record, 'size')
