import random


class Encoding:
    """Abstract base class for DNA encodings"""

    encoding_map = NotImplemented
    decoding_map = NotImplemented
    bits_per_base = NotImplemented

    @classmethod
    def encode(cls, bytes_encoded):
        """encode a DNA bytes string

        :param int bytes_encoded: encoded DNA string
        :return bytes: decoded bytes DNA sequence
        """
        raise NotImplementedError

    def decode(self, integer_encoded):
        """decode an encoded string

        :param int integer_encoded: encoded DNA string
        :return bytes: decoded bytes DNA sequence
        """
        raise NotImplementedError

    def gc_content(self, integer_encoded):
        raise NotImplementedError

    @staticmethod
    def hamming_distance(a, b):
        raise NotImplementedError


class TwoBit(Encoding):
    """Encode a DNA sequence using a 2-bit encoding.

    Two-bit encoding uses 0 for an encoded nucleotide. As such, it cannot distinguish between
    the end of sequence and trailing A nucleotides, and thus decoding these strings requires
    knowledge of their length. As such, it is only appropriate for encoding fixed sequence
    lengths

    In addition, in order to encode in 2-bit, N-nucleotides must be randomized to one of A, C,
    G, and T.

    :param int sequence_length: number of nucleotides that are being encoded
    """

    def __init__(self, sequence_length):
        self.sequence_length = sequence_length

    class TwoBitEncodingMap:
        """
        custom class that mimics a read-only dict, but generates random nucleotides for ambiguous
        nucleotides e.g. N
        """
        map_ = {ord('A'): 0, ord('C'): 1, ord('T'): 2, ord('G'): 3,
                ord('a'): 0, ord('c'): 1, ord('t'): 2, ord('g'): 3}

        iupac_ambiguous = {ord(c) for c in 'MRWSYKVHDBNmrwsykvhdbn'}

        def __getitem__(self, byte):
            try:
                return self.map_[byte]
            except KeyError:
                if byte not in self.iupac_ambiguous:
                    raise KeyError('%s is not a valid IUPAC nucleotide code' % chr(byte))
                return random.randint(0, 3)

    encoding_map = TwoBitEncodingMap()
    decoding_map = {0: b'A', 1: b'C', 2: b'T', 3: b'G'}
    bits_per_base = 2

    @classmethod
    def encode(cls, bytes_encoded):
        """
        encodes a DNA bytes string using 2 bits per base. N (and other ambiguous nucleotides)
        are randomized

        :param bytes bytes_encoded: a bytes string containing A, C, G, T nucleotides.
        :return int: encoded sequence
        """
        encoded = 0
        for character in bytes_encoded:
            encoded <<= 2
            encoded += cls.encoding_map[character]
        return encoded

    def decode(self, integer_encoded):
        """decodes an 2-bit encoded string

        :param int integer_encoded: encoded nucleotide sequence
        :return bytes: decoded sequence
        """
        decoded = b''
        for _ in range(self.sequence_length):
            decoded = self.decoding_map[integer_encoded & 3] + decoded
            integer_encoded >>= 2
        return decoded

    def gc_content(self, integer_encoded):
        """return the gc content of an 2-bit encoded bytes object

        :param int integer_encoded: encoded nucleotide string
        """
        i = 0
        for _ in range(self.sequence_length):
            i += integer_encoded & 1
            integer_encoded >>= 2
        return i

    @staticmethod
    def hamming_distance(a, b):
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

    def __init__(self, *args, **kwargs):
        """
        args and kwargs are not used, but allows ThreeBit to be initialized the same way as TwoBit,
        despite not requiring a sequence length parameter.
        """
        pass

    class ThreeBitEncodingMap:

        # C: 1, A: 2, G: 3, T: 4, N: 6;  # note, not using 0
        map_ = {ord('C'): 1, ord('A'): 2, ord('G'): 3, ord('T'): 4, ord('N'): 6,
                ord('c'): 1, ord('a'): 2, ord('g'): 3, ord('t'): 4, ord('n'): 6}

        def __getitem__(self, byte):
            try:
                return self.map_[byte]
            except KeyError:
                return 6  # any non-standard nucleotide gets "N"

    encoding_map = ThreeBitEncodingMap()
    decoding_map = {1: b'C', 2: b'A', 3: b'G', 4: b'T', 6: b'N'}
    bits_per_base = 3

    @classmethod
    def encode(cls, bytes_encoded):
        """
        encode a DNA bytes string using 3 bits per base.

        :param bytes bytes_encoded: a bytes string containing A, C, G, T, and N nucleotides.
        :return int: encoded string
        """
        encoded = 0
        for character in bytes_encoded:
            encoded <<= 3
            encoded += cls.encoding_map[character]
        return encoded

    @classmethod
    def decode(cls, integer_encoded):
        """decodes an 3-bit encoded string

        :param int integer_encoded: encoded nucleotide sequence
        :return bytes: decoded sequence
        """
        decoded = b''
        while integer_encoded:
            decoded = cls.decoding_map[integer_encoded & 7] + decoded
            integer_encoded >>= 3
        return decoded

    @classmethod
    def gc_content(cls, integer_encoded):
        """return the gc content of an 3-bit encoded bytes object

        :param int integer_encoded: encoded nucleotide string
        """
        i = 0
        while integer_encoded:
            i += integer_encoded & 1
            integer_encoded >>= 3
        return i

    @staticmethod
    def hamming_distance(a, b):
        difference = a ^ b
        d_hamming = 0
        while difference:
            if difference & 7:
                d_hamming += 1
            difference >>= 3
        return d_hamming
