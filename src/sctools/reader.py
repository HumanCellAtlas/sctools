import os
import gzip
import bz2
from copy import copy
from collections.abc import Iterable
from functools import partial
from typing import Callable


def infer_open(file_: str, mode: str) -> Callable:
    """
    Helper function to infer the correct inferred_openhook, for file_ ignoring extensions

    :param str file_: the file to open
    :param str mode: options: ['r', 'rb'] the intended open mode
    :return Callable: open function with mode pre-set through functools.partial
    """
    with open(file_, 'rb') as f:
        data: bytes = f.read(3)

        # gz and bzip treat 'r' = bytes, 'rt' = string
        if data[:2] == b'\x1f\x8b':  # gzip magic number
            inferred_openhook: Callable = gzip.open
            inferred_mode: str = 'rt' if mode == 'r' else mode

        elif data == b'BZh':  # bz2 magic number
            inferred_openhook: Callable = bz2.open
            inferred_mode: str = 'rt' if mode == 'r' else mode

        else:
            inferred_openhook: Callable = open
            inferred_mode: str = mode

    return partial(inferred_openhook, mode=inferred_mode)


class Reader:
    """
    Basic reader object that seamlessly loops over multiple input files

    Can be subclassed to create readers for specific file types (fastq, gtf, etc.)
    """

    def __init__(self, files='-', mode='r', header_comment_char=None):
        """
        :param list|str files: file or list of files to be read. Defaults to sys.stdin
        :param mode: (Default 'r') returns string objects. Change to 'rb' to
          return bytes objects.
        """

        if isinstance(files, str):
            self._files = [files]
        elif isinstance(files, Iterable):  # test items of iterable
            files = list(files)
            if all(isinstance(f, str) for f in files):
                self._files = files
            else:
                raise TypeError('all passed files must be type str')
        else:
            raise TypeError('files must be a string filename or a list of such names.')

        # set open mode:
        if mode not in {'r', 'rb'}:
            raise ValueError('mode must be one of r, rb')
        self._mode = mode

        if isinstance(header_comment_char, str) and mode == 'rb':
            self._header_comment_char = header_comment_char.encode()
        else:
            self._header_comment_char = header_comment_char

    @property
    def filenames(self):
        return self._files

    def __len__(self):
        """
        return the length of the Reader object.

        Note that this function requires reading the complete file, and should typically not be
        used with sys.stdin, as it will consume the input.
        """
        return sum(1 for _ in self)

    def __iter__(self):
        for file_ in self._files:

            f = infer_open(file_, self._mode)(file_)

            # iterate over the file, dropping header lines if requested
            try:
                file_iterator = iter(f)
                if self._header_comment_char is not None:
                    first_record = next(file_iterator)
                    while first_record.startswith(self._header_comment_char):
                        first_record = next(file_iterator)

                    yield first_record  # avoid loss of first non-comment line

                for record in file_iterator:  # now, run to exhaustion
                    yield record
            finally:  # clean up
                f.close()

    @property
    def size(self):
        """return the collective size of all files being read in bytes"""
        return sum(os.stat(f).st_size for f in self._files)

    def select_record_indices(self, indices):
        """iterate over provided indices only, skipping other records.

        :param set indices:
        :return Iterator:
        """
        indices = copy(indices)  # passed indices is a reference, need own copy to modify
        for idx, record in enumerate(self):
            if idx in indices:
                yield record
                indices.remove(idx)

                # stopping condition
                if not indices:
                    break


def zip_readers(*readers, indices=None):
    """zip together multiple reader objects, yielding records simultaneously.

    :param [Reader] readers:
    :param set indices: set of indices to iterate over

    :return Iterator: iterator over tuples of records, one from each passed Reader object.
    """
    if indices:
        iterators = zip(*(r.select_record_indices(indices) for r in readers))
    else:
        iterators = zip(*readers)
    for record_tuple in iterators:
        yield record_tuple
