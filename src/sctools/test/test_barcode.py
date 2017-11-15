import os
from .. import barcode, encodings
import numpy as np
import pytest

data_dir = os.path.split(__file__)[0] + '/data/'


@pytest.fixture
def barcode_set():
    return barcode.Barcodes.from_whitelist(
        data_dir + '1K-august-2016.txt', barcode_length=16)


@pytest.fixture(scope='module', params=['r', 'rb'])
def short_barcode_set_from_iterable(request):
    with open(data_dir + '1K-august-2016.txt', request.param) as f:
        barcodes = [l.strip() for l in f.readlines()[:50]]
    if isinstance(barcodes[0], bytes):
        return barcode.Barcodes.from_iterable_bytes(barcodes, barcode_length=16)
    else:
        return barcode.Barcodes.from_iterable_strings(barcodes, barcode_length=16)


@pytest.fixture(scope='module')
def short_barcode_set_from_encoded():
    return barcode.Barcodes.from_iterable_encoded([0, 1, 2, 3, 4, 5, 6, 7], barcode_length=2)


def test_iterable_produces_correct_barcodes(short_barcode_set_from_encoded):
    tbe = encodings.TwoBit(2)
    decoded = [tbe.decode(b) for b in short_barcode_set_from_encoded]
    print(decoded)
    assert decoded == [b'AA', b'AC', b'AT', b'AG', b'CA', b'CC', b'CT', b'CG']


def test_reads_barcodes_from_file(barcode_set):
    assert len(barcode_set) == 1000  # number of barcodes in file.


def test_base_frequency_sums_are_all_equal_to_barcode_set_length(barcode_set):
    bf = barcode_set.base_frequency()
    assert isinstance(bf, np.ndarray)
    assert np.array_equal(bf.sum(axis=1), np.ones(16) * len(barcode_set))


def test_barcode_diversity_is_in_range(barcode_set):
    bd = barcode_set.effective_diversity()
    assert np.all(bd >= 0)
    assert np.all(bd <= 1)


def test_summarize_hamming_distances_gives_reasonable_results(short_barcode_set_from_iterable):

    hamming_summary = short_barcode_set_from_iterable.summarize_hamming_distances()

    # we know 10x barcodes have at least this much distance
    assert hamming_summary['minimum'] >= 2
    # no barcode can have more hamming distance than length
    assert all(v <= 16 for v in hamming_summary.values())
