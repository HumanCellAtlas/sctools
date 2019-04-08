import os

import numpy as np
import pysam
import pytest

from .. import barcode, encodings, platform, consts

data_dir = os.path.split(__file__)[0] + '/data/'


# TEST BARCODES


@pytest.fixture
def barcode_set():
    return barcode.Barcodes.from_whitelist(
        data_dir + '1k-august-2016.txt', barcode_length=16
    )


@pytest.fixture(scope='module', params=['r', 'rb'])
def short_barcode_set_from_iterable(request):
    with open(data_dir + '1k-august-2016.txt', request.param) as f:
        barcodes = [l.strip() for l in f.readlines()[:50]]
    if isinstance(barcodes[0], bytes):
        return barcode.Barcodes.from_iterable_bytes(barcodes, barcode_length=16)
    else:
        return barcode.Barcodes.from_iterable_strings(barcodes, barcode_length=16)


@pytest.fixture(scope='module')
def short_barcode_set_from_encoded():
    return barcode.Barcodes.from_iterable_encoded(
        [0, 1, 2, 3, 4, 5, 6, 7], barcode_length=2
    )


def test_iterable_produces_correct_barcodes(short_barcode_set_from_encoded):
    tbe = encodings.TwoBit(2)
    decoded = [tbe.decode(b) for b in short_barcode_set_from_encoded]
    print(decoded)
    assert decoded == [b'AA', b'AC', b'AT', b'AG', b'CA', b'CC', b'CT', b'CG']


def test_reads_barcodes_from_file(barcode_set):
    assert len(barcode_set) == 1001  # number of barcodes in file.


def test_base_frequency_sums_are_all_equal_to_barcode_set_length(barcode_set):
    bf = barcode_set.base_frequency()
    assert isinstance(bf, np.ndarray)
    assert np.array_equal(bf.sum(axis=1), np.ones(16) * len(barcode_set))


def test_barcode_diversity_is_in_range(barcode_set):
    bd = barcode_set.effective_diversity()
    assert np.all(bd >= 0)
    assert np.all(bd <= 1)


def test_summarize_hamming_distances_gives_reasonable_results(
    short_barcode_set_from_iterable
):

    hamming_summary = short_barcode_set_from_iterable.summarize_hamming_distances()

    # we know 10x barcodes have at least this much distance
    assert hamming_summary['minimum'] >= 2
    # no barcode can have more hamming distance than length
    assert all(v <= 16 for v in hamming_summary.values())


# TEST HashErrorsToCorrectBarcodes


@pytest.fixture(scope='module')
def trivial_whitelist():
    barcode_iterable = ['A' * 8]
    error_mapping = barcode.ErrorsToCorrectBarcodesMap._prepare_single_base_error_hash_table(
        barcode_iterable
    )
    return barcode.ErrorsToCorrectBarcodesMap(error_mapping)


@pytest.fixture(scope='module')
def truncated_whitelist_from_10x():
    # note that this whitelist contains 1 non-10x barcode to ensure the presence of a matching
    # target in the test data.
    error_mapping = barcode.ErrorsToCorrectBarcodesMap.single_hamming_errors_from_whitelist(
        data_dir + '1k-august-2016.txt'
    )
    return error_mapping


def test_incorrect_input_raises_errors(trivial_whitelist):
    with pytest.raises(TypeError):
        barcode.ErrorsToCorrectBarcodesMap('not_a_mapping')
    with pytest.raises(TypeError):
        barcode.ErrorsToCorrectBarcodesMap({'not_a_mapping'})
    with pytest.raises(TypeError):
        barcode.ErrorsToCorrectBarcodesMap(['not_a_mapping', 'sldkf'])
    assert isinstance(trivial_whitelist, barcode.ErrorsToCorrectBarcodesMap)


def test_correct_barcode_finds_and_corrects_1_base_errors(trivial_whitelist):
    assert trivial_whitelist.get_corrected_barcode('TAAAAAAA') == 'AAAAAAAA'
    assert trivial_whitelist.get_corrected_barcode('AAAACAAA') == 'AAAAAAAA'
    assert trivial_whitelist.get_corrected_barcode('AAAGAAAA') == 'AAAAAAAA'
    assert trivial_whitelist.get_corrected_barcode('AAAAAAAA') == 'AAAAAAAA'


def test_correct_barcode_raises_keyerror_when_barcode_not_correct_length(
    trivial_whitelist
):
    with pytest.raises(KeyError):
        trivial_whitelist.get_corrected_barcode('AAA')
    with pytest.raises(KeyError):
        trivial_whitelist.get_corrected_barcode('AAAAAAAAA')
    with pytest.raises(KeyError):
        trivial_whitelist.get_corrected_barcode('AAAAAAAAAA')


def test_correct_barcode_raises_keyerror_when_barcode_has_more_than_one_error(
    trivial_whitelist
):
    with pytest.raises(KeyError):
        trivial_whitelist.get_corrected_barcode('AAAAAATT')
    with pytest.raises(KeyError):
        trivial_whitelist.get_corrected_barcode('TTAAAAAA')


@pytest.fixture(scope='module')
def tagged_bamfile():
    outbam = data_dir + 'bam_with_tags_test.bam'
    args = [
        '--r1',
        data_dir + 'test_r1.fastq',
        '--i1',
        data_dir + 'test_i7.fastq',
        '--u2',
        data_dir + 'test.bam',
        '--output-bamfile',
        outbam,
    ]
    platform.TenXV2.attach_barcodes(args)
    yield outbam
    os.remove(outbam)


def test_correct_bam_produces_cb_tags(tagged_bamfile, truncated_whitelist_from_10x):
    outbam = data_dir + 'bam_with_cb_tags.bam'
    truncated_whitelist_from_10x.correct_bam(tagged_bamfile, outbam)
    success = False
    with pysam.AlignmentFile(outbam, 'rb') as f:
        for record in f:
            try:
                success = record.get_tag(consts.CELL_BARCODE_TAG_KEY)
            except KeyError:
                continue
    assert success
    os.remove(outbam)
