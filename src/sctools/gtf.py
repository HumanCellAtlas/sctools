"""
GTF Records and Iterators
=========================

.. currentmodule:: sctools

This module defines a GTF record class and a Reader class to iterate over GTF-format files

Classes
-------
Record      Data class that exposes GTF record fields by name
Reader      GTF file reader that yields GTF Records

References
----------
https://useast.ensembl.org/info/website/upload/gff.html
"""

import string
from typing import List, Dict, Generator, Iterable

from . import reader


class Record:
    """Data class for storing and interacting with GTF records

    Subclassed to produce exon, transcript, and gene-specific record types.
    A gtf record has 8 fixed fields which are followed by optional fields separated by ;\t, which
    are stored by this class in the attributes field and accessible by get_attribute. Fixed fields
    are accessible by name.

    Parameters
    ----------
    record : str
        an unparsed GTF record

    Attributes
    ----------
    seqname : str
        The name of the sequence (often chromosome) this record is found on.
    chromosome : str
        Synonym for seqname.
    source : str
        The group responsible for generating this annotation.
    feature : str
        The type of record (e.g. gene, exon, ...).
    start : str
        The start position of this feature relative to the beginning of seqname.
    end : str
        The end position of this feature relative to the beginning of seqname....
    score : str
        The annotation score. Rarely used.
    strand : {'+', '-'}
        The strand of seqname that this annotation is found on
    frame : {'0', '1', '2'}
        '0' indicates that the first base of the feature is the first base of a codon,
        '1' that the second base is the first base of a codon, and so on
    size : int
        the number of nucleotides spanned by this feature

    Methods
    -------
    get_attribute(key: str)
        attempt to retrieve a variable field with name equal to `key`
    set_attribute(key: str, value: str)
        set variable field `key` equal to `value`. Overwrites `key` if already present.

    """

    __slots__ = ['_fields', '_attributes']

    _del_letters: str = string.ascii_letters
    _del_non_letters: str = ''.join(
        set(string.printable).difference(string.ascii_letters))

    def __init__(self, record: str):
        fields: List[str] = record.strip(';\n').split('\t')

        self._fields: List[str] = fields[:8]

        self._attributes: Dict[str, str] = {
            key: value.strip('"') for (key, value) in
            [field.strip().split() for field in fields[8].split(';')]
        }

    def __repr__(self):
        return '<Record: %s>' % self.__str__()

    def __bytes__(self):
        return self.__str__().encode()

    def __str__(self):
        return '\t'.join(self._fields) + self._format_attribute() + '\n'

    def __hash__(self) -> int:
        return hash(self.__str__())

    def _format_attribute(self):
        return ' '.join('%s "%s";' % (k, v) for k, v in self._attributes.items())

    @property
    def seqname(self) -> str:
        return self._fields[0]

    @property
    def chromosome(self) -> str:
        return self._fields[0]  # synonym for seqname

    @property
    def source(self) -> str:
        return self._fields[1]

    @property
    def feature(self) -> str:
        return self._fields[2]

    @property
    def start(self) -> int:
        return int(self._fields[3])

    @property
    def end(self) -> int:
        return int(self._fields[4])

    @property
    def score(self) -> str:
        return self._fields[5]

    @property
    def strand(self) -> str:
        return self._fields[6]

    @property
    def frame(self) -> str:
        return self._fields[7]

    @property
    def size(self) -> int:
        size = self.end - self.start
        if size < 0:
            raise ValueError('invalid record: negative size %d (start > end)' % size)
        else:
            return size

    def get_attribute(self, key) -> str:
        """access an item from the attribute field of a GTF file.

        Parameters
        ----------
        key : str
            Item to retrieve

        Returns
        -------
        value : str
            Contents of variable attribute `key`

        Raises
        ------
        KeyError
            if there is no variable attribute `key` associated with this record

        """
        return self._attributes.get(key)

    def set_attribute(self, key, value) -> None:
        """Set variable attribute `key` equal to `value`

        If attribute `key` is already set for this record, its contents are overwritten by `value`

        Parameters
        ----------
        key : str
            attribute name
        value : str
            attribute content

        """
        self._attributes[key] = value

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __ne__(self, other):
        return not self.__eq__(other)


class Reader(reader.Reader):
    """GTF file iterator

    Parameters
    ----------
    files : Union[str, List], optional
        File(s) to read. If '-', read sys.stdin (default = '-')
    mode : {'r', 'rb'}, optional
        Open mode. If 'r', read strings. If 'rb', read bytes (default = 'r').
    header_comment_char : str, optional
        lines beginning with this character are skipped (default = '#')

    Methods
    -------
    filter(retain_types: Iterable[str])
        Iterate over a gtf file, only yielding records in `retain_types`.
    __iter__()
        iterate over gtf records in file, yielding `Record` objects

    See Also
    --------
    sctools.reader.Reader

    """

    def __init__(self, files='-', mode='r', header_comment_char='#'):
        super().__init__(files, mode, header_comment_char)  # has different default args from super

    def __iter__(self):
        for line in super().__iter__():
            yield Record(line)

    def filter(self, retain_types: Iterable[str]) -> Generator:
        """Iterate over a gtf file, returning only record whose feature type is in retain_types.

        Features are stored in GTF field 2.

        Parameters
        ----------
        retain_types : Iterable[str]
            Record feature types to retain.

        Yields
        ------
        gtf_record : Record
            gtf `Record` object

        """
        retain_types = set(retain_types)
        for record in self:
            if record.feature in retain_types:
                yield record
