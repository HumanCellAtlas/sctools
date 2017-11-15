from collections import Iterable
import string
from . import reader


class Record:
    """
    Simple namespace object that makes the fields of a GTF record available. Subclassed
    to create records specific to exons, transcripts, and genes
    """

    __slots__ = ['_fields', '_attributes']

    _del_letters = string.ascii_letters
    _del_non_letters = ''.join(
        set(string.printable).difference(string.ascii_letters))

    def __init__(self, record):
        """
        :param str record: input record from file
        """
        # gtf has 8 fixed fields which are followed by optional fields separated by ;\t, which are
        # stored by this class in the attributes field.
        fields = record.strip(';\n').split('\t')
        self._fields = fields[:8]
        self._attributes = {
            key: value.strip('"') for (key, value) in
            [field.split() for field in fields[8].split('; ')]
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
    def seqname(self):
        return self._fields[0]

    @property
    def chromosome(self):
        return self._fields[0]  # synonym for seqname

    @property
    def source(self):
        return self._fields[1]

    @property
    def feature(self):
        return self._fields[2]

    @property
    def start(self):
        return int(self._fields[3])

    @property
    def end(self):
        return int(self._fields[4])

    @property
    def score(self):
        return self._fields[5]

    @property
    def strand(self):
        return self._fields[6]

    @property
    def frame(self):
        return self._fields[7]

    @property
    def size(self):
        size = self.end - self.start
        if size < 0:
            raise ValueError('invalid record: negative size %d (start > end)' % size)
        else:
            return size

    def get_attribute(self, key):
        """access an item from the attribute field of a GTF file.

        :param key: item to access
        :return str: value of item
        """
        return self._attributes.get(key)

    def set_attribute(self, key, value):
        """set a gtf attribute

        :param str key: attribute name
        :param str value: attribute value
        """
        self._attributes[key] = value

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __ne__(self, other):
        return not self.__eq__(other)


class Reader(reader.Reader):
    """
    SubClass of reader.Reader, returns an Reader with several specialized iterator
    methods.

    :method __iter__: Iterator over all non-header records in gtf; yields Record objects.
    :method iter_genes: Iterator over all genes in gtf; yields Gene objects.
    """

    def __init__(self, files='-', mode='r', header_comment_char='#'):
        """
        :param list|str files: (default sys.stdin) file or list of files to be read.
        :param mode: (Default 'r') returns string objects. Change to 'rb' to
          return bytes.
        :param str|bytes header_comment_char: (default #) character that marks headers, which are
          ignored
        """
        super().__init__(files, mode, header_comment_char)  # has different default args from super

    def __iter__(self):
        for line in super().__iter__():
            yield Record(line)

    def filter(self, retain_types):
        """iterate over a gtf file, returning only record whose feature type is in retain_types.

        :param Iterable retain_types: a set of record feature types to retain

        :return Iterator:
        """
        retain_types = set(retain_types)
        for record in self:
            if record.feature in retain_types:
                yield record
