import pytest
from .. import encodings
from itertools import combinations


@pytest.fixture(scope="module")
def sequence():
    return b"ACGTTTGAGATGAGATATAGANNNN"


@pytest.fixture(scope="module")
def encoder_2bit(sequence):
    length = len(sequence)
    return encodings.TwoBit(length)


@pytest.fixture(scope="module")
def encoder_3bit():
    return encodings.ThreeBit()


@pytest.fixture(scope="module", params=[encodings.TwoBit, encodings.ThreeBit])
def encoder(request):
    return request.param


def test_two_bit_encode_decode_produces_same_string_except_for_N(
    sequence, encoder_2bit
):
    encoded = encoder_2bit.encode(sequence)
    decoded = encoder_2bit.decode(encoded)
    assert sequence[:4] == decoded[:4]  # last 4 are N, which get randomized


def test_three_bit_encode_decode_produces_same_string(sequence, encoder_3bit):
    encoded = encoder_3bit.encode(sequence)
    decoded = encoder_3bit.decode(encoded)
    assert sequence == decoded


def test_two_bit_encoder_gets_correct_gc_content(encoder_2bit):
    sequence_no_n = b"AGCGCGAT"
    gc_content = sequence_no_n.count(b"C") + sequence_no_n.count(b"G")
    encoded = encoder_2bit.encode(sequence_no_n)
    assert encoder_2bit.gc_content(encoded) == gc_content


def test_three_bit_encoder_gets_correct_gc_content(sequence, encoder_3bit):
    encoded = encoder_3bit.encode(sequence)
    assert encoder_3bit.gc_content(encoded) == sequence.count(b"C") + sequence.count(
        b"G"
    )


def test_two_bit_throws_errors_when_asked_to_encode_unknown_nucleotide(encoder_2bit):
    with pytest.raises(KeyError):
        encoder_2bit.encode(b"ACGTP")  # P is not a valid code


def test_three_bit_encodes_unknown_nucleotides_as_N(encoder_3bit):
    encoded = encoder_3bit.encode(b"ACGTP")  # P is not a valid code
    decoded = encoder_3bit.decode(encoded)
    assert decoded == b"ACGTN"


@pytest.fixture
def simple_barcodes():
    """simple barcode set with min_hamming = 1, max_hamming = 2"""
    return [b"ACGT", b"ACGG", b"ACGA", b"ACGC", b"TCGT", b"CCGT", b"GCGT"]


@pytest.fixture
def simple_hamming_distances(simple_barcodes):
    simple_hamming_distances = []
    for a, b in combinations(simple_barcodes, 2):
        d_hamming = 0
        for i, j in zip(a, b):
            if i != j:
                d_hamming += 1
        simple_hamming_distances.append(d_hamming)
    return simple_hamming_distances


def test_encoded_hamming_distance_is_accurate(
    simple_hamming_distances, simple_barcodes, encoder
):
    # encode simple barcodes
    tbe = encoder(4)
    encoded = [tbe.encode(b) for b in simple_barcodes]
    encoded_hamming_distances = []

    # use hamming distance function
    for a, b in combinations(encoded, 2):
        encoded_hamming_distances.append(tbe.hamming_distance(a, b))

    # verify they are the same as the simple function used in this file
    assert simple_hamming_distances == encoded_hamming_distances
