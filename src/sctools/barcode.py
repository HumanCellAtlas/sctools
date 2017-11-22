import itertools
from collections import Counter, Mapping
import numpy as np
from .encodings import TwoBit
from .stats import base4_entropy
import pysam


class Barcodes:

    def __init__(self, barcodes, barcode_length):
        """
        container for a set of barcodes, encoded in 2bit representation for fast operations. Can
        optionally be constructed from an iterable where barcodes can be present multiple times. In
        these cases, barcodes are analyzed based on their observed frequencies.

        :param Mapping barcodes: dictionary
        :param int barcode_length: the length of all barcodes in the set. Different-length barcodes
          are not supported
        """
        if not isinstance(barcodes, Mapping):
            raise TypeError('barcode set must be a dict-like object mapping barcodes to counts')
        self._mapping = barcodes
        if not isinstance(barcode_length, int) and barcode_length > 0:
            raise ValueError('barcode length must be a positive integer')
        self._barcode_length = barcode_length

    def __contains__(self, item):
        return item in self._mapping

    def __iter__(self):
        return iter(self._mapping)

    def __len__(self):
        return len(self._mapping)

    def __getitem__(self, item):
        return self._mapping[item]

    def summarize_hamming_distances(self):
        """returns descriptive statistics on hamming distances between pairs of barcodes"""
        distances = []
        for a, b in itertools.combinations(self, 2):
            distances.append(TwoBit.hamming_distance(a, b))
        keys = ('minimum', '25th percentile', 'median', '75th percentile', 'maximum', 'average')
        values = list(np.percentile(distances, [0, 25, 50, 75, 100])) + [np.mean(distances)]
        return dict(zip(keys, values))

    def base_frequency(self, weighted=False):
        """return the frequency of each base at each position in the barcode set

        Note: weighting is currently not supported, and must be set to False or will raise
         NotImplementedError

        :param bool weighted: (Default False) if True, each barcode is counted once for each time it
          was observed
        :return np.array: barcode_length x 4 2d array
        """
        base_counts_by_position = np.zeros((self._barcode_length, 4), dtype=np.uint64)
        keys = np.fromiter(self._mapping.keys(), dtype=np.uint64)

        for i in reversed(range(self._barcode_length)):
            binary_base_representations, counts = np.unique(keys & 3, return_counts=True)
            if weighted:  # todo weighted not working, values multiplication does not work
                raise NotImplementedError
            else:
                base_counts_by_position[i, binary_base_representations] = counts
            keys >>= 2

        return base_counts_by_position

    def effective_diversity(self, weighted=False):
        """returns the effective base diversity of the barcode set by position.

        maximum diversity for each position is 1, and represents a perfect split of 25% per base at
        a given position.

        :param bool weighted: (Default False) if True, each barcode is counted once for each time it
          was observed

        :return np.array[float]: array of size barcode_length containing floats in [0, 1]
        """
        return base4_entropy(self.base_frequency(weighted=weighted))

    @classmethod
    def from_whitelist(cls, file_, barcode_length):
        """Creates a barcode set from a whitelist file

        :param str file_: location of the whitelist file. Should be formatted one barcode per line.
          barcodes should be encoded in plain text (UTF-8, ASCII), not bit-encoded.
          each barcode will be assigned a count of 1.
        :param int barcode_length: length of the barcodes in the file.

        :return BarcodeSet:
        """
        tbe = TwoBit(barcode_length)
        with open(file_, 'rb') as f:
            return cls(Counter(tbe.encode(barcode[:-1]) for barcode in f), barcode_length)

    @classmethod
    def from_iterable_encoded(cls, iterable, barcode_length):
        """construct an ObservedBarcodeSet from an iterable of encoded barcodes"""
        return cls(Counter(iterable), barcode_length=barcode_length)

    @classmethod
    def from_iterable_strings(cls, iterable, barcode_length):
        """construct an ObservedBarcodeSet from an iterable of string barcodes"""
        tbe = TwoBit(barcode_length)
        return cls(Counter(tbe.encode(b.encode()) for b in iterable), barcode_length=barcode_length)

    @classmethod
    def from_iterable_bytes(cls, iterable, barcode_length):
        """construct an ObservedBarcodeSet from an iterable of bytes barcodes"""
        tbe = TwoBit(barcode_length)
        return cls(Counter(tbe.encode(b) for b in iterable), barcode_length=barcode_length)


class ErrorsToCorrectBarcodesMap:

    def __init__(self, errors_to_barcodes):
        """use a hash map of erroneous_barcode -> true barcode to correct errors in cell barcodes

        :param dict errors_to_barcodes: mapping of errors to the barcodes that could generate them.
        """
        if not isinstance(errors_to_barcodes, Mapping):
            raise TypeError('errors_to_barcodes must be a mapping of erroneous barcodes to correct '
                            'barcodes, not %a' % type(errors_to_barcodes))
        self._map = errors_to_barcodes

    def get_corrected_barcode(self, barcode):
        return self._map[barcode]

    @staticmethod
    def _prepare_single_base_error_hash_table(barcodes):
        """

        :param Iterable barcodes: iterable of string barcodes
        :return dict: mapping between erroneous barcodes with single-base mutations and the barcode
          they were generated from
        """
        error_map = {}
        for barcode in barcodes:

            # include correct barcode
            error_map[barcode] = barcode

            # include all single-base errors
            for i, nucleotide in enumerate(barcode):
                errors = set('ACGTN')
                errors.discard(nucleotide)
                for e in errors:
                    error_map[barcode[:i] + e + barcode[i + 1:]] = barcode
        return error_map

    @classmethod
    def single_hamming_errors_from_whitelist(cls, whitelist_file):
        """factory method to generate instance of class from a file containing "correct" barcodes

        :param str whitelist_file: text file containing barcode per line
        :return ErrorsToCorrectBarcodesMap: instance of cls.
        """
        with open(whitelist_file, 'r') as f:
            return cls(cls._prepare_single_base_error_hash_table((line[:-1] for line in f)))

    def correct_bam(self, bam_file, output_bam_file):
        """correct barcodes in a (potentially unaligned) bamfile, given a whitelist

        :param str bam_file: BAM format file in same order as the fastq files
        :param str output_bam_file:  BAM format file containing cell, umi, and sample tags.
        """
        with pysam.AlignmentFile(bam_file, 'rb') as fin, \
                pysam.AlignmentFile(output_bam_file, 'wb', template=fin) as fout:
            for alignment in fin:
                try:
                    # pysam tags mimic tags in SAM specification; Z = string
                    tag = self.get_corrected_barcode(alignment.get_tag('CR'))
                except KeyError:  # pass through the uncorrected barcode.
                    tag = alignment.get_tag('CR')
                alignment.set_tag(tag='CB', value=tag, value_type='Z')
                fout.write(alignment)
