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

import logging
import string
from typing import List, Dict, Generator, Iterable, Union

from . import reader

_logger = logging.getLogger(__name__)


class GTFRecord:
    """Data class for storing and interacting with GTF records

    Subclassed to produce exon, transcript, and gene-specific record types.
    A GTF record has 8 fixed fields which are followed by optional fields separated by ;\t, which
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
        set(string.printable).difference(string.ascii_letters)
    )

    def __init__(self, record: str):
        fields: List[str] = record.strip(';\n').split('\t')

        self._fields: List[str] = fields[:8]

        self._attributes: Dict[str, str] = {}
        for field in fields[8].split(';'):
            try:
                key, _, value = field.strip().partition(' ')
                self._attributes[key] = value.strip('"')
            except Exception:
                raise RuntimeError(
                    f'Error parsing field "{field}" of GTF record "{record}"'
                )

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
            raise ValueError(f'Invalid record: negative size {size} (start > end)')
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
        Iterate over a GTF file, only yielding records in `retain_types`.
    __iter__()
        iterate over GTF records in file, yielding `Record` objects

    See Also
    --------
    sctools.reader.Reader

    """

    def __init__(self, files='-', mode='r', header_comment_char='#'):
        super().__init__(
            files, mode, header_comment_char
        )  # has different default args from super

    def __iter__(self):
        for line in super().__iter__():
            yield GTFRecord(line)

    def filter(self, retain_types: Iterable[str]) -> Generator:
        """Iterate over a GTF file, returning only record whose feature type is in retain_types.

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


# todo this lenient behavior is deemed to change in the future (warning -> exception)
def _resolve_multiple_gene_names(gene_name: str):
    _logger.warning(
        f'Multiple entries encountered for "{gene_name}". Please validate the input GTF file(s). '
        f'Skipping the record for now; in the future, this will be considered as a '
        f'malformed GTF file.'
    )


def extract_gene_names(
    files: Union[str, List[str]] = '-', mode: str = 'r', header_comment_char: str = '#'
) -> Dict[str, int]:
    """Extract gene names from GTF file(s) and returns a map from gene names to their corresponding
    occurrence orders in the given file(s).

    Parameters
    ----------
    files : Union[str, List], optional
        File(s) to read. If '-', read sys.stdin (default = '-')
    mode : {'r', 'rb'}, optional
        Open mode. If 'r', read strings. If 'rb', read bytes (default = 'r').
    header_comment_char : str, optional
        lines beginning with this character are skipped (default = '#')

    Returns
    -------
    Dict[str, int]
        A map from gene names to their linear index
    """
    gene_name_to_index: Dict[str, int] = dict()
    gene_index = 0
    for record in Reader(files, mode, header_comment_char).filter(
        retain_types=['gene']
    ):
        gene_name = record.get_attribute('gene_name')
        if gene_name is None:
            raise ValueError(
                f'Malformed GTF file detected. Record is of type gene but does not have a '
                f'"gene_name" field: {record}'
            )
        if gene_name in gene_name_to_index:
            _resolve_multiple_gene_names(gene_name)
            continue
        gene_name_to_index[gene_name] = gene_index
        gene_index += 1
    return gene_name_to_index

def extract_extended_gene_names(
    files: Union[str, List[str]] = '-', mode: str = 'r', header_comment_char: str = '#'
) -> Dict[str, tuple]:
    """Extract extended gene names from GTF file(s) and returns a map from gene names to their corresponding
    occurrence locations the given file(s).

    Parameters
    ----------
    files : Union[str, List], optional
        File(s) to read. If '-', read sys.stdin (default = '-')
    mode : {'r', 'rb'}, optional
        Open mode. If 'r', read strings. If 'rb', read bytes (default = 'r').
    header_comment_char : str, optional
        lines beginning with this character are skipped (default = '#')

    Returns
    -------
    List[tuple(int,int), key]
        A map from gene names to their start and end tuples
    """
    gene_name_to_start_end: Dict[str, List[tuple]] = dict()
    for record in Reader(files, mode, header_comment_char).filter(
        retain_types=['gene']
    ):
        gene_name = record.get_attribute('gene_name')
        if gene_name is None:
            raise ValueError(
                f'Malformed GTF file detected. Record is of type gene but does not have a '
                f'"gene_name" field: {record}'
            )
        if gene_name in gene_name_to_start_end:
            _resolve_multiple_gene_names(gene_name)
            continue

        if not record.chromosome in gene_name_to_start_end: 
           gene_name_to_start_end[record.chromosome] = dict()
        gene_name_to_start_end[record.chromosome][gene_name] = (record.start, record.end)

    gene_locations = dict()
    for chromosome in gene_name_to_start_end: 
        gene_locations[chromosome] = [ (locs, key) for key, locs in  gene_name_to_start_end[chromosome].items() ] 
        gene_locations[chromosome].sort(key =lambda x: x[0])

    return gene_locations
