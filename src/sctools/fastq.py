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
mock_10x_i1(input_fastq, index_barcode, index_quality)
    creates an i1 file that matches an r1 or r2 `input_fastq` file

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

from . import reader
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
            raise TypeError('fastq name must be bytes')
        elif not value.startswith(b'@'):
            raise ValueError('fastq name must start with @')
        else:
            self._name = value

    @property
    def sequence(self) -> AnyStr:
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        """fastq nucleotide sequence"""
        if not isinstance(value, (bytes, str)):
            raise TypeError('fastq sequence must be str or bytes')
        else:
            self._sequence = value

    @property
    def name2(self) -> AnyStr:
        return self._name2

    @name2.setter
    def name2(self, value):
        """second fastq record name field (rarely used)"""
        if not isinstance(value, (bytes, str)):
            raise TypeError('fastq name2 must be str or bytes')
        else:
            self._name2 = value

    @property
    def quality(self) -> AnyStr:
        return self._quality

    @quality.setter
    def quality(self, value):
        """fastq record base call quality scores"""
        if not isinstance(value, (bytes, str)):
            raise TypeError('fastq quality must be str or bytes')
        else:
            self._quality = value

    def __bytes__(self):
        return b''.join((self.name, self.sequence, self.name2, self.quality))

    def __str__(self):
        return b''.join((self.name, self.sequence, self.name2, self.quality)).decode()

    def __repr__(self):
        return "Name: %s\nSequence: %s\nName2: %s\nQuality: %s\n" % (
            self.name, self.sequence, self.name2, self.quality)

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
        Iterable of 4 bytes strings that comprise a fastq record

    Attributes
    ----------
    name : str
        fastq record name
    sequence : str
        fastq nucleotide sequence
    name2 : str
        second fastq record name field (rarely used)
    quality : str
        base call quality for each nucleotide in sequence

    Methods
    -------
    average_quality()
        The average quality of the fastq record

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
        """fastq record name"""
        if not isinstance(value, (bytes, str)):
            raise TypeError('fastq name must be str or bytes')
        if not value.startswith('@'):
            raise ValueError('fastq name must start with @')
        else:
            self._name = value

    def average_quality(self) -> float:
        """return the average quality of this record"""
        b = self.quality[:-1].encode()
        return sum(c for c in b) / len(b) - 33  # -33 due to solexa/illumina phred conversion


class Reader(reader.Reader):
    """Fastq Reader that defines some special methods for reading and summarizing fastq data.

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
        """Iterate over a fastq file, returning records

        Yields
        ------
        fastq_record : Tuple[str]
            tuple of length 4 containing the name, sequence, name2, and quality for a fastq record

        """
        record_type = StrRecord if self._mode == 'r' else Record
        for record in self._record_grouper(super().__iter__()):
            yield record_type(record)


# namedtuple that defines the start and end position of a barcode sequence and provides the name
# for both a quality and sequence tag
EmbeddedBarcode = namedtuple('Tag', ['start', 'end', 'sequence_tag', 'quality_tag'])


def extract_barcode(record, embedded_barcode) -> Tuple[Tuple[str, str, str], Tuple[str, str, str]]:
    """Extracts barcodes from a fastq record at positions defined by an EmbeddedBarcode object.

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
    seq = record.sequence[embedded_barcode.start:embedded_barcode.end]
    qual = record.quality[embedded_barcode.start:embedded_barcode.end]
    return (embedded_barcode.sequence_tag, seq, 'Z'), (embedded_barcode.quality_tag, qual, 'Z')


# todo the reader subclasses need better docs
class EmbeddedBarcodeGenerator(Reader):
    """Generate barcodes from a fastq file(s) from positions defined by EmbeddedBarcode(s)

    Extracted barcode objects are produced in a form that is consumable by pysam's bam and sam
    set_tag methods.

    Parameters
    ----------
    embedded_barcodes : Iterable[EmbeddedBarcode]
        tag objects defining start and end of the sequence containing the tag, and the tag
        identifiers for sequence and quality tags
    fastq_files : str | List, optional
        fastq file or files to be read. (default = sys.stdin)
    mode : {'r', 'rb'}, optional
        open mode for fastq files. If 'r', return string. If 'rb', return bytes (default = 'r')

    """

    def __init__(self, fastq_files, embedded_barcodes, *args, **kwargs):
        super().__init__(files=fastq_files, *args, **kwargs)
        self.embedded_barcodes = embedded_barcodes

    def __iter__(self):
        """iterates over barcodes extracted from fastq"""
        for record in super().__iter__():  # iterates records; we extract barcodes.
            barcodes = []
            for barcode in self.embedded_barcodes:
                barcodes.extend(extract_barcode(record, barcode))
            yield barcodes


# todo the reader subclasses need better docs
class BarcodeGeneratorWithCorrectedCellBarcodes(Reader):
    """Generate barcodes from fastq file(s) from positions defined by EmbeddedBarcode(s)

    Extracted barcode objects are produced in a form that is consumable by pysam's bam and sam
    set_tag methods. In this class, one EmbeddedBarcode must be defined as an
    `embedded_cell_barcode`, which is checked against a whitelist and error corrected during
    generation

    Parameters
    ----------
    fastq_files : str | List, optional
        fastq file or files to be read. (default = sys.stdin)
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
            other_embedded_barcodes: Iterable[EmbeddedBarcode]=tuple(),
            *args, **kwargs):

        super().__init__(files=fastq_files, *args, **kwargs)
        if isinstance(other_embedded_barcodes, (list, tuple)):
            self.embedded_barcodes = other_embedded_barcodes
        else:
            raise TypeError('if passed, other_embedded_barcodes must be a list or tuple')

        self._error_mapping = ErrorsToCorrectBarcodesMap.single_hamming_errors_from_whitelist(
            whitelist)
        self.embedded_cell_barcode = embedded_cell_barcode

    def __iter__(self):
        """iterates over barcodes extracted from fastq"""
        for record in super().__iter__():  # iterates records; we extract barcodes.
            barcodes = []

            barcodes.extend(self.extract_cell_barcode(record, self.embedded_cell_barcode))
            for barcode in self.embedded_barcodes:
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
        seq_tag, qual_tag = extract_barcode(record, cb)
        try:
            corrected_cb = self._error_mapping.get_corrected_barcode(seq_tag[1])
            return seq_tag, qual_tag, ('CB', corrected_cb, 'Z')
        except KeyError:
            return seq_tag, qual_tag

def mock_10x_i1(input_fastq: str, index_barcode: str, index_quality: str) -> None:
    """Generate a synthetic i1 file that matches the read names in input_fastq

    The I1 file will be output in the same directory as the input file with a correctly
    matched filename for a 10x experiment

    Parameters
    ----------
    input_fastq: str
        filename of r1 or r2 fastq file to match read names of
    index_barcode: str
        the sample index that will be used for all I1 reads
    index_quality: str
        the sample quality that will be used for all I1 reads

    """

    if not sample_barcode.endswith('\n'):
        sample_barcode += '\n'
    if not sample_quality.endswith('\n'):
        sample_quality += '\n'

    output_filename = re.sub('_[r|R][1|2]', '_I1_', input_fastq)
    rd = sctools.fastq.Reader(input_fastq)
    with open(output_filename, 'w') as f:
        for record in rd:
            record.sequence = sample_barcode
            record.quality = sample_quality
            f.write(str(record))

