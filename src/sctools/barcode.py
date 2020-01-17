"""
Nucleotide Barcode Manipulation Tools
=====================================

.. currentmodule:: sctools

This module contains tools to characterize oligonucleotide barcodes and a simple hamming-base
error-correction approach which corrects barcodes within a specified distance of a "whitelist" of
expected barcodes.

Classes
-------
Barcodes                        Class to characterize a set of barcodes
ErrorsToCorrectBarcodesMap      Class to carry out error correction routines

"""

import itertools
from collections import Counter
from typing import Mapping, Iterator, List, Tuple, Iterable

import numpy as np
import pysam

from . import consts
from .encodings import TwoBit
from .stats import base4_entropy


class Barcodes:
    """Container for a set of nucleotide barcodes.

    Contained barcodes are encoded in 2bit representation for fast operations. Instances of this
    class can optionally be constructed from an iterable where barcodes can be present multiple
    times. In these cases, barcodes are analyzed based on their observed frequencies.

    Parameters
    ----------
    barcodes: Mapping[str, int]
        dictionary-like mapping barcodes to the number of times they were observed
    barcode_length: int
        the length of all barcodes in the set. Different-length barcodes are not supported.

    See Also
    --------
    sctools.encodings.TwoBit

    """

    def __init__(self, barcodes: Mapping[str, int], barcode_length: int):
        if not isinstance(barcodes, Mapping):
            raise TypeError(
                'The argument "barcodes" must be a dict-like object mapping barcodes to counts'
            )
        self._mapping: Mapping[str, int] = barcodes

        if not isinstance(barcode_length, int) and barcode_length > 0:
            raise ValueError('The argument "barcode_length" must be a positive integer')
        self._barcode_length: int = barcode_length

    def __contains__(self, item) -> bool:
        return item in self._mapping

    def __iter__(self) -> Iterator[str]:
        return iter(self._mapping)

    def __len__(self) -> int:
        return len(self._mapping)

    def __getitem__(self, item) -> int:
        return self._mapping[item]

    def summarize_hamming_distances(self) -> Mapping[str, float]:
        """Returns descriptive statistics on hamming distances between pairs of barcodes.

        Returns
        -------
        descriptive_statistics : Mapping[str, float]
            minimum, 25th percentile, median, 75th percentile, maximum, and average hamming
            distance between all pairs of barcodes

        References
        ----------
        https://en.wikipedia.org/wiki/Hamming_distance

        """
        distances: List = []

        for a, b in itertools.combinations(self, 2):
            distances.append(TwoBit.hamming_distance(a, b))

        keys: Tuple = (
            "minimum",
            "25th percentile",
            "median",
            "75th percentile",
            "maximum",
            "average",
        )
        values: List = list(np.percentile(distances, [0, 25, 50, 75, 100]))
        values.append(np.mean(distances))

        return dict(zip(keys, values))

    def base_frequency(self, weighted=False) -> np.ndarray:
        """return the frequency of each base at each position in the barcode set

        Notes
        -----
        weighting is currently not supported, and must be set to False or base_frequency will raise
        NotImplementedError  # todo fix

        Parameters
        ----------
        weighted: bool, optional
            if True, each barcode is counted once for each time it was observed (default = False)

        Returns
        -------
        frequencies : np.array
            barcode_length x 4 2d numpy array

        Raises
        ------
        NotImplementedError
            if weighted is True

        """
        base_counts_by_position: np.ndarray = np.zeros(
            (self._barcode_length, 4), dtype=np.uint64
        )

        keys: np.ndarray = np.fromiter(self._mapping.keys(), dtype=np.uint64)

        for i in reversed(range(self._barcode_length)):
            binary_base_representations, counts = np.unique(
                keys & 3, return_counts=True
            )
            if weighted:
                raise NotImplementedError
            else:
                base_counts_by_position[i, binary_base_representations] = counts

            # finished with this nulceotide, move two bits forward to the next one
            keys >>= 2

        return base_counts_by_position

    def effective_diversity(self, weighted=False) -> np.ndarray:
        """Returns the effective base diversity of the barcode set by position.

        maximum diversity for each position is 1, and represents a perfect split of 25% per base at
        a given position.

        Parameters
        ----------
        weighted : bool, optional
            if True, each barcode is counted once for each time it was observed (default = False)

        Returns
        -------
        effective_diversity : np.array[float]
            1-d array of size barcode_length containing floats in [0, 1]

        """
        return base4_entropy(self.base_frequency(weighted=weighted))

    @classmethod
    def from_whitelist(cls, file_: str, barcode_length: int):
        """Creates a barcode set from a whitelist file.

        Parameters
        ----------
        file_ : str
            location of the whitelist file. Should be formatted one barcode per line. Barcodes
            should be encoded in plain text (UTF-8, ASCII), not bit-encoded. Each barcode will be
            assigned a count of 1.
        barcode_length : int
            Length of the barcodes in the file.

        Returns
        -------
        barcodes : Barcodes
            class object containing barcodes from a whitelist file

        """
        tbe = TwoBit(barcode_length)
        with open(file_, "rb") as f:
            return cls(
                Counter(tbe.encode(barcode[:-1]) for barcode in f), barcode_length
            )

    @classmethod
    def from_iterable_encoded(cls, iterable: Iterable[int], barcode_length: int):
        """Construct an ObservedBarcodeSet from an iterable of encoded barcodes.

        Parameters
        ----------
        iterable : Iterable[int]
            iterable of barcodes encoded in TwoBit representation
        barcode_length : int
            the length of the barcodes in `iterable`

        Returns
        -------
        barcodes : Barcodes
            class object containing barcodes from a whitelist file
        """
        return cls(Counter(iterable), barcode_length=barcode_length)

    @classmethod
    def from_iterable_strings(cls, iterable: Iterable[str], barcode_length: int):
        """Construct an ObservedBarcodeSet from an iterable of string barcodes.

        Parameters
        ----------
        iterable : Iterable[str]
            iterable of barcodes encoded in TwoBit representation
        barcode_length : int
            the length of the barcodes in `iterable`

        Returns
        -------
        barcodes : Barcodes
            class object containing barcodes from a whitelist file
        """
        tbe: TwoBit = TwoBit(barcode_length)
        return cls(
            Counter(tbe.encode(b.encode()) for b in iterable),
            barcode_length=barcode_length,
        )

    @classmethod
    def from_iterable_bytes(cls, iterable: Iterable[bytes], barcode_length: int):
        """Construct an ObservedBarcodeSet from an iterable of bytes barcodes.

        Parameters
        ----------
        iterable : Iterable[bytes]
            iterable of barcodes in bytes representation
        barcode_length : int
            the length of the barcodes in `iterable`

        Returns
        -------
        barcodes : Barcodes
            class object containing barcodes from a whitelist file
        """
        tbe: TwoBit = TwoBit(barcode_length)
        return cls(
            Counter(tbe.encode(b) for b in iterable), barcode_length=barcode_length
        )


class ErrorsToCorrectBarcodesMap:
    """Correct any barcode that is within one hamming distance of a whitelisted barcode

    Parameters
    ----------
    errors_to_barcodes : Mapping[str, str]
        dict-like mapping 1-base errors to the whitelist barcode that they could be generated from

    Methods
    -------
    get_corrected_barcode(barcode: str)
        Return a barcode if it is whitelist, or the corrected version if within edit distance 1
    correct_bam(bam_file: str, output_bam_file: str)
        correct barcodes in a bam file, given a whitelist

    References
    ----------
    https://en.wikipedia.org/wiki/Hamming_distance

    """

    def __init__(self, errors_to_barcodes: Mapping[str, str]):
        if not isinstance(errors_to_barcodes, Mapping):
            raise TypeError(
                f'The argument "errors_to_barcodes" must be a mapping of erroneous barcodes to correct '
                f"barcodes, not {type(errors_to_barcodes)}"
            )
        self._map = errors_to_barcodes

    def get_corrected_barcode(self, barcode: str) -> str:
        """Return a barcode if it is whitelist, or the corrected version if within edit distance 1

        Parameters
        ----------
        barcode : str
            the barcode to return the corrected version of. If the barcode is in the whitelist,
            the input barcode is returned unchanged.

        Returns
        -------
        corrected_barcode : str
            corrected version of the barcode

        Raises
        ------
        KeyError
            if the passed barcode is not within 1 hamming distance of any whitelist barcode

        References
        ----------
        https://en.wikipedia.org/wiki/Hamming_distance

        """
        return self._map[barcode]

    @staticmethod
    def _prepare_single_base_error_hash_table(
        barcodes: Iterable[str],
    ) -> Mapping[str, str]:
        """Generate a map of correct barcodes and single base error codes to whitelist barcodes

        Parameters
        ----------
        barcodes : Iterable[str]
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
                errors = set("ACGTN")
                errors.discard(nucleotide)
                for e in errors:
                    error_map[barcode[:i] + e + barcode[i + 1 :]] = barcode
        return error_map

    @classmethod
    def single_hamming_errors_from_whitelist(cls, whitelist_file: str):
        """Factory method to generate instance of class from a file containing "correct" barcodes.

        Parameters
        ----------
        whitelist_file : str
            Text file containing barcode per line.

        Returns
        -------
        errors_to_barcodes_map : ErrorsToCorrectBarcodesMap
            instance of cls, built from whitelist

        """
        with open(whitelist_file, "r") as f:
            return cls(
                cls._prepare_single_base_error_hash_table((line[:-1] for line in f))
            )

    def correct_bam(self, bam_file: str, output_bam_file: str) -> None:
        """Correct barcodes in a (potentially unaligned) bamfile, given a whitelist.

        Parameters
        ----------
        bam_file : str
            BAM format file in same order as the fastq files
        output_bam_file : str
            BAM format file containing cell, umi, and sample tags.

        """
        with pysam.AlignmentFile(bam_file, "rb") as fin, pysam.AlignmentFile(
            output_bam_file, "wb", template=fin
        ) as fout:
            for alignment in fin:
                try:
                    tag = self.get_corrected_barcode(alignment.get_tag("CR"))
                except KeyError:  # pass through the uncorrected barcode.
                    tag = alignment.get_tag(consts.RAW_CELL_BARCODE_TAG_KEY)
                alignment.set_tag(
                    tag=consts.CELL_BARCODE_TAG_KEY, value=tag, value_type="Z"
                )
                fout.write(alignment)
