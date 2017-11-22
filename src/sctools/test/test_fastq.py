from itertools import product
from functools import partial
import string
import os
import pytest
from .. import fastq
from ..reader import zip_readers


# set some useful globals for testing
data_dir = os.path.split(__file__)[0] + '/data/'
_i7_files = [data_dir + f for f in ('test_i7.fastq', 'test_i7.fastq.gz', 'test_i7.fastq.bz2')]
_files = [data_dir + f for f in ('test_i7.fastq', 'test_r1.fastq', 'test_r2.fastq')]
_gz_files = [data_dir + f for f in ('test_i7.fastq.gz', 'test_r1.fastq.gz', 'test_r2.fastq.gz')]
_bz2_files = [data_dir + f for f in ('test_i7.fastq.bz2', 'test_r1.fastq.bz2', 'test_r2.fastq.bz2')]

_modes = ('r', 'rb')
_files_and_modes = list(product(_i7_files, _modes))
_multifiles_and_modes = list(product((_files, _gz_files, _bz2_files), _modes))
_map_encoder = {'r': str, 'rb': partial(bytes, encoding='utf-8')}


# TEST READER

@pytest.fixture(scope='module', params=_files_and_modes)
def i7_files_compressions_and_modes(request):
    """generates different compression types and modes for testing"""
    return request.param[0], request.param[1]


@pytest.fixture(scope='module', params=_multifiles_and_modes)
def reader_all_compressions(request):
    """generates open fastq reader files for each compression and read mode"""
    return fastq.Reader(request.param[0], request.param[1])


@pytest.fixture(scope='module')
def bytes_fastq_record():
    return [
        b'@name\n',
        b'ACTACAAT\n',
        b'+\n',
        b'%%%%AAAA\n']


@pytest.fixture(scope='module')
def string_fastq_record():
    return [
        '@name\n',
        'ACTACAAT\n',
        '+\n',
        '%%%%AAAA\n']


def test_reader_stores_filenames():
    names = ['notreal', 'fake']
    rd = fastq.Reader(files=names)
    assert rd.filenames == names


def test_reader_reads_first_record(reader_all_compressions):
    for record in reader_all_compressions:
        assert isinstance(record, fastq.Record)
        expected_result = 'NCACAATG\n' if isinstance(record.sequence, str) else b'NCACAATG\n'
        assert record.sequence == expected_result
        break  # just first record


def test_reader_skips_header_character_raises_value_error(i7_files_compressions_and_modes):
    """
    test should skip the first name line, shifting each record up 1. As a result, the
     first sequence should be found in the name field
    """
    filename, mode = i7_files_compressions_and_modes
    rd = fastq.Reader(
        filename,
        mode=mode,
        header_comment_char='@')

    with pytest.raises(ValueError):
        next(iter(rd))


def test_reader_reads_correct_number_of_records_across_multiple_files(reader_all_compressions):
    assert len(reader_all_compressions) == 300  # 3 files


def test_mixed_filetype_read_gets_correct_record_number():
    rd = fastq.Reader(
        [_gz_files[0], _bz2_files[0]],
        mode='r',
        header_comment_char='#')

    assert len(rd) == 200


def test_non_string_filename_raises_typeerror():
    with pytest.raises(TypeError):
        _ = fastq.Reader(10, 'r')


def test_non_string_filename_in_iterable_raises_typeerror():
    with pytest.raises(TypeError):
        _ = fastq.Reader(('works', 10), 'r')


def test_invalid_open_mode_raises_valueerror():
    with pytest.raises(ValueError):
        _ = fastq.Reader('works', 'not_acceptable_open_mode')


def test_fastq_returns_correct_filesize_for_single_and_multiple_files():
    rd = fastq.Reader(
        _i7_files[0],
        mode='r',  # mode irrelevant
        header_comment_char='#')
    assert rd.size == 7774

    rd = fastq.Reader(
        _i7_files,
        mode='r',  # mode irrelevant
        header_comment_char='#')
    assert rd.size == 7774 + 853 + 802  # three file sizes


def test_reader_properly_subsets_based_on_indices():
    rd = fastq.Reader(_i7_files[0], mode='r')
    indices = {0, 5, 10, 12}
    n_records = sum(1 for _ in rd.select_record_indices(indices))
    assert n_records == len(indices)


def test_zipping_readers_generates_expected_output():
    rd1 = fastq.Reader(_files[0], 'r')
    rd2 = fastq.Reader(_files[0], 'r')
    for r1, r2 in zip_readers(rd1, rd2):
        assert isinstance(r1, fastq.Record)
        assert isinstance(r2, fastq.Record)
        expected_result = 'NCACAATG\n'
        assert r1.sequence == r2.sequence == expected_result
        break  # just first record


def test_zipping_readers_with_indices_generates_expected_output():
    rd1 = fastq.Reader(_files[0], 'r')
    rd2 = fastq.Reader(_files[0], 'r')
    indices = {0, 1, 2, 3}
    for r1, r2 in zip_readers(rd1, rd2, indices=indices):
        assert isinstance(r1, fastq.Record)
        assert isinstance(r2, fastq.Record)
        expected_result = 'NCACAATG\n'
        assert r1.sequence == r2.sequence == expected_result
        break  # just first record


def test_printing_bytes_record_generates_valid_fastq_record(bytes_fastq_record):
    record = fastq.Record(bytes_fastq_record)
    assert str(record) == b''.join(bytes_fastq_record).decode()
    assert bytes(record) == b''.join(bytes_fastq_record)


def test_bytes_fastq_record_quality_score_parsing(bytes_fastq_record):
    record = fastq.Record(bytes_fastq_record)
    assert record.average_quality() == 18


def test_printing_string_record_generates_valid_fastq_record(string_fastq_record):
    record = fastq.StrRecord(string_fastq_record)
    assert str(record) == ''.join(string_fastq_record)
    assert bytes(record) == ''.join(string_fastq_record).encode()


def test_string_fastq_record_quality_score_parsing(string_fastq_record):
    record = fastq.StrRecord(string_fastq_record)
    assert record.average_quality() == 18


# TEST RECORD

def test_fields_populate_properly(reader_all_compressions):
    encoder = _map_encoder[reader_all_compressions._mode]
    name_prefix = encoder('@')
    alphabet = set(encoder('ACGTN'))
    name2_string = encoder('+\n')
    ascii_chars = set(i for i in encoder(string.printable))
    for record in reader_all_compressions:
        assert record.name.startswith(name_prefix)
        assert all(i in alphabet for i in record.sequence.strip())
        assert record.name2 == name2_string
        assert all(i in ascii_chars for i in record.quality.strip())


# TEST BarcodeGeneratorWithCorrectedCellbarcodes

@pytest.fixture(scope='function')
def embedded_barcode_generator():
    cell_barcode = fastq.EmbeddedBarcode(start=0, end=16, quality_tag='CY', sequence_tag='CR')
    molecule_barcode = fastq.EmbeddedBarcode(start=16, end=26, quality_tag='UY', sequence_tag='UR')
    return fastq.EmbeddedBarcodeGenerator(data_dir + 'test_r1.fastq.gz',
                                          [cell_barcode, molecule_barcode])


@pytest.fixture(scope='function')
def barcode_generator_with_corrected_cell_barcodes():
    cell_barcode = fastq.EmbeddedBarcode(start=0, end=16, quality_tag='CY', sequence_tag='CR')
    molecule_barcode = fastq.EmbeddedBarcode(start=16, end=26, quality_tag='UY', sequence_tag='UR')
    return fastq.BarcodeGeneratorWithCorrectedCellBarcodes(
        data_dir + 'test_r1.fastq.gz', cell_barcode, data_dir + '1k-august-2016.txt',
        [molecule_barcode])


def test_embedded_barcode_generator_produces_outputs_of_expected_size(embedded_barcode_generator):
    for cell_seq, cell_qual, umi_seq, umi_qual in embedded_barcode_generator:

        # correct values
        correct_cell_barcode_length = 16
        correct_umi_length = 10

        # note that all barcodes are strings and therefore should get 'Z' values

        # test cell tags
        assert cell_seq[0] == 'CR'
        assert len(cell_seq[1]) == correct_cell_barcode_length
        assert all(v in 'ACGTN' for v in cell_seq[1])
        assert cell_seq[2] == 'Z'
        assert cell_qual[0] == 'CY'
        assert len(cell_qual[1]) == correct_cell_barcode_length
        assert all(v in string.printable for v in cell_qual[1])
        assert cell_seq[2] == 'Z'

        # test umi tags
        assert umi_seq[0] == 'UR'
        assert len(umi_seq[1]) == correct_umi_length
        assert all(v in 'ACGTN' for v in umi_seq[1])
        assert umi_seq[2] == 'Z'
        assert umi_qual[0] == 'UY'
        assert len(umi_qual[1]) == correct_umi_length
        assert all(v in string.printable for v in umi_qual[1])
        assert umi_seq[2] == 'Z'

        break  # just the first tag is fine


def test_corrects_barcodes(barcode_generator_with_corrected_cell_barcodes):
    success = False
    for barcode_sets in barcode_generator_with_corrected_cell_barcodes:
        for barcode_set in barcode_sets:
            if barcode_set[0] == 'CB':
                success = True
                break
    assert success
