"""
Compressed Barcode Encoding Methods
===================================

.. currentmodule:: sctools

This module defines several classes to encode DNA sequences in memory-efficient forms, using 2 bits
to encode bases of a 4-letter DNA alphabet (ACGT) or 3 bits to encode a 5-letter DNA alphabet
that includes the ambiguous call often included by Illumina base calling software (ACGTN). The
classes also contain several methods useful for efficient querying and manipulation of the encoded
sequence.

Classes
-------
Encoding                            Encoder base class
ThreeBit                            Three bit DNA encoder / decoder
TwoBit                              Two bit DNA encoder / decoder

"""

import random
from typing import Mapping, AnyStr, Set


class Encoding:
    """

    Attributes
    ----------
    encoding_map : TwoBitEncodingMap
        Class that mimics a Mapping[bytes, str] where bytes must be a single byte encoded character
        (encoder)
    decoding_map : Mapping[int, bytes]
        Dictionary that maps integers to bytes human-readable representations (decoder)
    bits_per_base : int
        number of bits used to encode each base

    Methods
    -------
    encode(bytes_encoded: bytes)
        encode a DNA string in a compressed representation
    decode(integer_encoded: int)
        decode a compressed DNA string into a human readable bytes format
    gc_content(integer_encoded: int)
        calculate the GC content of an encoded DNA string
    hamming_distance(a: int, b: int)
        calculate the hamming distance between two encoded DNA strings

    """

    encoding_map: Mapping[AnyStr, int] = NotImplemented
    decoding_map: Mapping[int, AnyStr] = NotImplemented
    bits_per_base: int = NotImplemented

    @classmethod
    def encode(cls, bytes_encoded: bytes) -> int:
        """Encode a DNA bytes string.

        Parameters
        ----------
        bytes_encoded : bytes
            bytes DNA string

        Returns
        -------
        encoded : int
            Encoded DNA sequence

        """
        raise NotImplementedError

    def decode(self, integer_encoded: int) -> bytes:
        """Decode a DNA bytes string.

        Parameters
        ----------
        integer_encoded : bytes
            Integer encoded DNA string

        Returns
        -------
        decoded : bytes
            Bytes decoded DNA sequence

        """
        raise NotImplementedError

    def gc_content(self, integer_encoded: int) -> int:
        """Return the number of G or C nucleotides in `integer_encoded`

        Parameters
        ----------
        integer_encoded : int
            Integer encoded DNA string

        Returns
        -------
        gc_content, int
            number of bases in `integer_encoded` input that are G or C.

        """
        raise NotImplementedError

    @staticmethod
    def hamming_distance(a, b) -> int:
        """Calculate the hamming distance between two DNA sequences

        The hamming distance counts the number of bases that are not the same nucleotide

        Parameters
        ----------
        a, b : int
            integer encoded


        Returns
        -------
        d : int
            hamming distance between a and b
        """
        raise NotImplementedError


class TwoBit(Encoding):
    """Encode a DNA sequence using a 2-bit encoding.

    Two-bit encoding uses 0 for an encoded nucleotide. As such, it cannot distinguish between
    the end of sequence and trailing A nucleotides, and thus decoding these strings requires
    knowledge of their length. Therefore, it is only appropriate for encoding fixed sequence
    lengths

    In addition, in order to encode in 2-bit, N-nucleotides must be randomized to one of A, C,
    G, and T.

    Parameters
    ----------
    sequence_length : int
        number of nucleotides that are being encoded

    """

    __doc__ += Encoding.__doc__

    def __init__(self, sequence_length: int):
        self.sequence_length: int = sequence_length

    class TwoBitEncodingMap:
        """Dict-like class that maps bytes to 2-bit integer representations

        Generates random nucleotides for ambiguous nucleotides e.g. N

        """

        map_ = {
            ord("A"): 0,
            ord("C"): 1,
            ord("T"): 2,
            ord("G"): 3,
            ord("a"): 0,
            ord("c"): 1,
            ord("t"): 2,
            ord("g"): 3,
        }

        iupac_ambiguous: Set[int] = {ord(c) for c in "MRWSYKVHDBNmrwsykvhdbn"}

        def __getitem__(self, byte: int) -> int:
            try:
                return self.map_[byte]
            except KeyError:
                if byte not in self.iupac_ambiguous:
                    raise KeyError(f"{chr(byte)} is not a valid IUPAC nucleotide code")
                return random.randint(0, 3)

    encoding_map: TwoBitEncodingMap = TwoBitEncodingMap()
    decoding_map: Mapping[int, bytes] = {0: b"A", 1: b"C", 2: b"T", 3: b"G"}
    bits_per_base: int = 2

    @classmethod
    def encode(cls, bytes_encoded: bytes) -> int:
        encoded = 0
        for character in bytes_encoded:
            encoded <<= 2
            encoded += cls.encoding_map[character]
        return encoded

    def decode(self, integer_encoded: int) -> bytes:
        decoded = b""
        for _ in range(self.sequence_length):
            decoded = self.decoding_map[integer_encoded & 3] + decoded
            integer_encoded >>= 2
        return decoded

    def gc_content(self, integer_encoded: int) -> int:
        i = 0
        for _ in range(self.sequence_length):
            i += integer_encoded & 1
            integer_encoded >>= 2
        return i

    @staticmethod
    def hamming_distance(a: int, b: int) -> int:
        difference = a ^ b
        d_hamming = 0
        while difference:
            if difference & 3:
                d_hamming += 1
            difference >>= 2
        return d_hamming


class ThreeBit(Encoding):
    """Encode a DNA sequence using a 3-bit encoding.

    Since no bases are encoded as 0, an empty triplet is interpreted as the end of the encoded
    string; Three-bit encoding can be used to encode and decode strings without knowledge of their
    length.

    """

    __doc__ += Encoding.__doc__

    def __init__(self, *args, **kwargs):
        """
        Notes
        -----
        args and kwargs are not used, but allow ThreeBit to be initialized the same way as TwoBit,
        despite not requiring a sequence length parameter.

        """
        pass

    class ThreeBitEncodingMap:
        """Dict-like class that maps bytes to 3-bit integer representations

        All IUPAC ambiguous codes are treated as "N"

        """

        # C: 1, A: 2, G: 3, T: 4, N: 6;  # note, not using 0
        map_ = {
            ord("C"): 1,
            ord("A"): 2,
            ord("G"): 3,
            ord("T"): 4,
            ord("N"): 6,
            ord("c"): 1,
            ord("a"): 2,
            ord("g"): 3,
            ord("t"): 4,
            ord("n"): 6,
        }

        def __getitem__(self, byte: int) -> int:
            try:
                return self.map_[byte]
            except KeyError:
                return 6  # any non-standard nucleotide gets "N"

    encoding_map: ThreeBitEncodingMap = ThreeBitEncodingMap()
    decoding_map: Mapping[int, bytes] = {1: b"C", 2: b"A", 3: b"G", 4: b"T", 6: b"N"}
    bits_per_base: int = 3

    @classmethod
    def encode(cls, bytes_encoded: bytes) -> int:
        encoded = 0
        for character in bytes_encoded:
            encoded <<= 3
            encoded += cls.encoding_map[character]
        return encoded

    @classmethod
    def decode(cls, integer_encoded: int) -> bytes:
        decoded = b""
        while integer_encoded:
            decoded = cls.decoding_map[integer_encoded & 7] + decoded
            integer_encoded >>= 3
        return decoded

    @classmethod
    def gc_content(cls, integer_encoded: int) -> int:
        i = 0
        while integer_encoded:
            i += integer_encoded & 1
            integer_encoded >>= 3
        return i

    @staticmethod
    def hamming_distance(a: int, b: int) -> int:
        difference = a ^ b
        d_hamming = 0
        while difference:
            if difference & 7:
                d_hamming += 1
            difference >>= 3
        return d_hamming
