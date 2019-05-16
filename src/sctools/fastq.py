"""
Efficient Fastq Iterators and Representations
=============================================

.. currentmodule:: sctools

This module implements classes for representing fastq records, reading and writing them, and
extracting parts of fastq sequence for transformation into bam format tags

Methods
-------
extract_barcode(record, embedded_barcode)
    extract a barcode, defined by `embedded_barcode` from `record`

Classes
-------
Record                                      Represents fastq records (input as bytes)
StrRecord                                   Represents fastq records (input as str)
Reader                                      Opens and iterates over fastq files
EmbeddedBarcodeGenerator                    Generates barcodes from a fastq file
BarcodeGeneratorWithCorrectedCellBarcodes   Generates (corrected) barcodes from a fastq file

References
----------
https://en.wikipedia.org/wiki/FASTQ_format

"""

from collections import namedtuple
from typing import Iterable, AnyStr, Iterator, Union, Tuple
import re

from . import reader, consts
from .barcode import ErrorsToCorrectBarcodesMap


# todo the inheritance pattern of this class is a bit confusing, particularly the str vs. bytes
# in the daughter classes
class Record:
    """Fastq Record.

    Parameters
    ----------
    record : Iterable[bytes]
        Iterable of 4 bytes strings that comprise a fastq record

    Attributes
    ----------
    name : bytes
        fastq record name
    sequence : bytes
        fastq nucleotide sequence
    name2 : bytes
        second fastq record name field (rarely used)
    quality : bytes
        base call quality for each nucleotide in sequence

    Methods
    -------
    average_quality()
        The average quality of the fastq record

    """

    __slots__ = ['_name', '_sequence', '_name2', '_quality']

    def __init__(self, record: Iterable[AnyStr]):
        # use the setter functions
        self.name, self.sequence, self.name2, self.quality = record

    @property
    def name(self) -> AnyStr:
        return self._name

    @name.setter
    def name(self, value):
        """fastq record name"""
        if not isinstance(value, (bytes, str)):
            raise TypeError('FASTQ name must be bytes')
        elif not value.startswith(b'@'):
            raise ValueError('FASTQ name must start with @')
        else:
            self._name = value

    @property
    def sequence(self) -> AnyStr:
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        """FASTQ nucleotide sequence"""
        if not isinstance(value, (bytes, str)):
            raise TypeError('FASTQ sequence must be str or bytes')
        else:
            self._sequence = value

    @property
    def name2(self) -> AnyStr:
        return self._name2

    @name2.setter
    def name2(self, value):
        """second FASTQ record name field (rarely used)"""
        if not isinstance(value, (bytes, str)):
            raise TypeError('FASTQ name2 must be str or bytes')
        else:
            self._name2 = value

    @property
    def quality(self) -> AnyStr:
        return self._quality

    @quality.setter
    def quality(self, value):
        """FASTQ record base call quality scores"""
        if not isinstance(value, (bytes, str)):
            raise TypeError('FASTQ quality must be str or bytes')
        else:
            self._quality = value

    def __bytes__(self):
        return b''.join((self.name, self.sequence, self.name2, self.quality))

    def __str__(self):
        return b''.join((self.name, self.sequence, self.name2, self.quality)).decode()

    def __repr__(self):
        return "Name: %s\nSequence: %s\nName2: %s\nQuality: %s\n" % (
            self.name,
            self.sequence,
            self.name2,
            self.quality,
        )

    def __len__(self):
        return len(self.sequence)

    def average_quality(self) -> float:
        """return the average quality of this record"""
        # -33 due to solexa/illumina phred conversion
        return sum(c for c in self.quality[:-1]) / (len(self.quality) - 1) - 33


class StrRecord(Record):
    """Fastq Record.

    Parameters
    ----------
    record : Iterable[str]
        Iterable of 4 bytes strings that comprise a FASTQ record

    Attributes
    ----------
    name : str
        FASTQ record name
    sequence : str
        FASTQ nucleotide sequence
    name2 : str
        second FASTQ record name field (rarely used)
    quality : str
        base call quality for each nucleotide in sequence

    Methods
    -------
    average_quality()
        The average quality of the FASTQ record

    """

    def __bytes__(self):
        return ''.join((self.name, self.sequence, self.name2, self.quality)).encode()

    def __str__(self):
        return ''.join((self.name, self.sequence, self.name2, self.quality))

    # todo is this method necessary?
    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, value):
        """FASTQ record name"""
        if not isinstance(value, (bytes, str)):
            raise TypeError('FASTQ name must be str or bytes')
        if not value.startswith('@'):
            raise ValueError('FASTQ name must start with @')
        else:
            self._name = value

    def average_quality(self) -> float:
        """return the average quality of this record"""
        b = self.quality[:-1].encode()
        return (
            sum(c for c in b) / len(b) - 33
        )  # -33 due to solexa/illumina phred conversion


class Reader(reader.Reader):
    """Fastq Reader that defines some special methods for reading and summarizing FASTQ data.

    Simple reader class that exposes an __iter__ and __len__  method

    Examples
    --------
    #todo add examples

    See Also
    --------
    sctools.reader.Reader

    References
    ----------
    https://en.wikipedia.org/wiki/FASTQ_format

    """

    @staticmethod
    def _record_grouper(iterable):
        """Groups contents of an iterator, yielding 4 objects at a time instead of one

        This is a somewhat complex python function. It creates 4 iterators on the same iterable;
        each moves the pointer to the position in the iterable forward when called, yielding 4
        objects at a time

        Returns
        -------
        grouped_iterator : Iterator[Str], Iterator[Str], Iterator[Str], Iterator[Str]

        """
        args = [iter(iterable)] * 4
        return zip(*args)

    def __iter__(self) -> Iterator[Tuple[str]]:
        """Iterate over a FASTQ file, returning records

        Yields
        ------
        fastq_record : Tuple[str]
            tuple of length 4 containing the name, sequence, name2, and quality for a FASTQ record

        """
        record_type = StrRecord if self._mode == 'r' else Record
        for record in self._record_grouper(super().__iter__()):
            yield record_type(record)


# namedtuple that defines the start and end position of a barcode sequence and provides the name
# for both a quality and sequence tag
class EmbeddedBarcode:
    def __init__(self, start, end, sequence_tag, quality_tag,
                 variable_length=False, molecular=True, minimum_length=8):
        self.start = start
        self.sequence_tag = sequence_tag
        self.quality_tag = quality_tag
        self.variable_length = variable_length
        self.end = end
        self.minimum_length = minimum_length


class VariableEmbeddedBarcode(EmbeddedBarcode) :
    def __init__(self, cell_sequence_tag, cell_quality_tag, molecule_sequence_tag, molecule_quality_tag,
                 minimum_length=8):
        self.cell_sequence_tag = cell_sequence_tag
        self.cell_quality_tag = cell_quality_tag
        self.molecule_sequence_tag = molecule_sequence_tag
        self.molecule_quality_tag = molecule_quality_tag
        self.minimum_length = minimum_length

#EmbeddedBarcode = namedtuple('Tag', ['start', 'end', 'sequence_tag', 'quality_tag'])


def extract_barcode(
    record, embedded_barcode
) -> Tuple[Tuple[str, str, str], Tuple[str, str, str]]:
    """Extracts barcodes from a FASTQ record at positions defined by an EmbeddedBarcode object.

    Parameters
    ----------
    record : FastqRecord
        Record to extract from
    embedded_barcode : EmbeddedBarcode
        Defines the barcode start and end positions and the tag name for the sequence and quality
        tags

    Returns
    -------
    sequence_tag : Tuple[str, str, 'Z']
        sequence tag identifier, sequence, SAM tag type ('Z' implies a string tag)
    quality_tag : Tuple[str, str, 'Z']
        quality tag identifier, quality, SAM tag type ('Z' implies a string tag)

    """
    seq = record.sequence[embedded_barcode.start : embedded_barcode.end]
    qual = record.quality[embedded_barcode.start : embedded_barcode.end]
    return (
        (embedded_barcode.sequence_tag, seq, 'Z'),
        (embedded_barcode.quality_tag, qual, 'Z'),
    )


def extract_variable_barcode(
    record, embedded_barcode
) -> Tuple[Tuple[str, str, str], Tuple[str, str, str]]:

    match = re.search(consts.SPLITTER_REGEX, record.sequence[embedded_barcode.minimum_length:])
    if match:
        cb_end_1 = match.start() + embedded_barcode.minimum_length + 1
    else:
        cb_end_1 = 8

    cb_start_2 = cb_end_1 + 22
    cb_end_2 = cb_start_2 + 8
    mb_start = cb_end_2
    mb_end = mb_start + 6

    cb_seq = record.sequence[0:cb_end_1] + record.sequence[cb_start_2:cb_end_2]
    cb_qual = record.quality[0:cb_end_1] + record.quality[cb_start_2:cb_end_2]
    mb_seq = record.sequence[mb_start:mb_end]
    mb_qual = record.quality[mb_start:mb_end]

    return(
        (embedded_barcode.cell_sequence_tag, cb_seq, 'Z'),
        (embedded_barcode.cell_quality_tag, cb_qual, 'Z'),
        (embedded_barcode.molecule_sequence_tag, mb_seq, 'Z'),
        (embedded_barcode.molecule_quality_tag, mb_qual, 'Z'),
    )


# todo the reader subclasses need better docs
class EmbeddedBarcodeGenerator(Reader):
    """Generate barcodes from a FASTQ file(s) from positions defined by EmbeddedBarcode(s)

    Extracted barcode objects are produced in a form that is consumable by pysam's bam and sam
    set_tag methods.

    Parameters
    ----------
    embedded_barcodes : Iterable[EmbeddedBarcode]
        tag objects defining start and end of the sequence containing the tag, and the tag
        identifiers for sequence and quality tags
    fastq_files : str | List, optional
        FASTQ file or files to be read. (default = sys.stdin)
    mode : {'r', 'rb'}, optional
        open mode for FASTQ files. If 'r', return string. If 'rb', return bytes (default = 'r')

    """

    def __init__(self, fastq_files, embedded_barcodes, extract_barcode_function=extract_barcode, *args, **kwargs):
        super().__init__(files=fastq_files, *args, **kwargs)
        self.embedded_barcodes = embedded_barcodes
        self.extract_barcode_function = extract_barcode_function

    def __iter__(self):
        """iterates over barcodes extracted from FASTQ"""
        for record in super().__iter__():  # iterates records; we extract barcodes.
            barcodes = []
            for barcode in self.embedded_barcodes:
                barcodes.extend(self.extract_barcode_function(record, barcode))
            yield barcodes


# todo the reader subclasses need better docs
class BarcodeGeneratorWithCorrectedCellBarcodes(Reader):
    """Generate barcodes from FASTQ file(s) from positions defined by EmbeddedBarcode(s)

    Extracted barcode objects are produced in a form that is consumable by pysam's bam and sam
    set_tag methods. In this class, one EmbeddedBarcode must be defined as an
    `embedded_cell_barcode`, which is checked against a whitelist and error corrected during
    generation

    Parameters
    ----------
    fastq_files : str | List, optional
        FASTQ file or files to be read. (default = sys.stdin)
    mode : {'r', 'rb'}, optional
        open mode for fastq files. If 'r', return string. If 'rb', return bytes (default = 'r')
    whitelist : str
        whitelist file containing "correct" cell barcodes for an experiment
    embedded_cell_barcodes : EmbeddedBarcode
        EmbeddedBarcode containing information about the position and names of cell barcode tags
    other_embedded_barcodes : Iterable[EmbeddedBarcode], optional
        tag objects defining start and end of the sequence containing the tag, and the tag
        identifiers for sequence and quality tags (default = None)

    Methods
    -------
    extract_cell_barcode(record: Record, cb: str)

    """

    def __init__(
        self,
        fastq_files: Union[str, Iterable[str]],
        embedded_cell_barcode: EmbeddedBarcode,
        whitelist: str,
        other_embedded_barcodes: Iterable[EmbeddedBarcode] = tuple(),
        is_variable=False,
        whitelist_2=None,
        *args,
        **kwargs
    ):

        super().__init__(files=fastq_files, *args, **kwargs)
        if isinstance(other_embedded_barcodes, (list, tuple)):
            self.embedded_barcodes = other_embedded_barcodes
        else:
            raise TypeError(
                'if passed, other_embedded_barcodes must be a list or tuple'
            )

        if is_variable:
            self._error_mapping = ErrorsToCorrectBarcodesMap.single_hamming_errors_from_inDrop_whitelist(
                whitelist, whitelist_2
            )
        else:
            self._error_mapping = ErrorsToCorrectBarcodesMap.single_hamming_errors_from_whitelist(
                whitelist
            )
        self.embedded_cell_barcode = embedded_cell_barcode
        self.is_variable = is_variable

    def __iter__(self):
        """iterates over barcodes extracted from fastq"""
        for record in super().__iter__():  # iterates records; we extract barcodes.
            barcodes = []

            barcodes.extend(
                self.extract_cell_barcode(record, self.embedded_cell_barcode)
            )
            for barcode in self.embedded_barcodes:
                if self.is_variable :
                    barcodes.extend(extract_variable_barcode(record, barcode))
                else :
                    barcodes.extend(extract_barcode(record, barcode))

            yield barcodes

    def extract_cell_barcode(self, record: Tuple[str], cb: EmbeddedBarcode):
        """Extract a cell barcode from a fastq record

        Parameters
        ----------
        record : Tuple[str]
            fastq record comprised of four strings: name, sequence, name2, and quality
        cb : EmbeddedBarcode
            defines the position and tag identifier for a call barcode

        Returns
        -------
        sequence_tag : Tuple[str, str, 'Z']
            raw sequence tag identifier, sequence, SAM tag type ('Z' implies a string tag)
        quality_tag : Tuple[str, str, 'Z']
            quality tag identifier, quality, SAM tag type ('Z' implies a string tag)
        corrected_tag : Optional[Tuple[str, str, 'Z']]
            Whitelist verified sequence tag. Only present if the raw sequence tag is in the
            whitelist or within 1 hamming distance of one of its barcodes

        """
        if self.is_variable :
            seq_tag, qual_tag, seq_tag_2, qual_tag_2 = extract_variable_barcode(record, cb)
            try:
                corrected_cb = self._error_mapping.get_corrected_barcode(seq_tag[1])
                return seq_tag, qual_tag, (consts.CELL_BARCODE_TAG_KEY, corrected_cb, 'Z'), seq_tag_2, qual_tag_2
            except KeyError:
                return seq_tag, qual_tag, seq_tag_2, qual_tag_2
        else :
            seq_tag, qual_tag = extract_barcode(record, cb)
            try:
                corrected_cb = self._error_mapping.get_corrected_barcode(seq_tag[1])
                return seq_tag, qual_tag, (consts.CELL_BARCODE_TAG_KEY, corrected_cb, 'Z')
            except KeyError:
                return seq_tag, qual_tag
