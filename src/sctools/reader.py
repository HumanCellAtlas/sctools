"""
Sequence File Iterators
=======================

.. currentmodule:: sctools

This module defines a general iterator and some helper functions for iterating over files
that contain sequencing data

Methods
-------
infer_open(file_: str, mode: str)
    helper function that determines the compression type of a file without relying on its extension
zip_readers(*readers, indices=None)
    helper function that iterates over one or more readers, optionally extracting only the records
    that correspond to indices

Classes
-------
Reader          Basic reader that loops over one or more input files.

See Also
--------
sctools.gtf.Reader
sctools.fastq.Reader

"""

import os
import gzip
import bz2
from copy import copy
from functools import partial
from typing import Callable, Iterable, Generator, Set, List


def infer_open(file_: str, mode: str) -> Callable:
    """Helper function to infer the correct compression type of an input file

    Identifies files that are .gz or .bz2 compressed without requiring file extensions

    Parameters
    ----------
    file_ : str
        the file to open
    mode : {'r', 'rb'}
        the mode to open the file in. 'r' returns strings, 'rb' returns bytes

    Returns
    -------
    open_function : Callable
        the correct open function for the file's compression with mode pre-set through functools
        partial

    """
    with open(file_, "rb") as f:
        data: bytes = f.read(3)

        # gz and bzip treat 'r' = bytes, 'rt' = string
        if data[:2] == b"\x1f\x8b":  # gzip magic number
            inferred_openhook: Callable = gzip.open
            inferred_mode: str = "rt" if mode == "r" else mode

        elif data == b"BZh":  # bz2 magic number
            inferred_openhook: Callable = bz2.open
            inferred_mode: str = "rt" if mode == "r" else mode

        else:
            inferred_openhook: Callable = open
            inferred_mode: str = mode

    return partial(inferred_openhook, mode=inferred_mode)


class Reader:
    """Basic reader object that seamlessly loops over multiple input files.

    Is subclassed to create readers for specific file types (e.g. fastq, gtf, etc.)

    Parameters
    ----------
    files : Union[str, List], optional
        The file(s) to read. If '-', read sys.stdin (default = '-')
    mode : {'r', 'rb'}, optional
        The open mode for files. If 'r', yield string data, if 'rb', yield bytes data
        (default = 'r').
    header_comment_char : str, optional
        If not None, skip lines beginning with this character (default = None).

    """

    def __init__(self, files="-", mode="r", header_comment_char=None):
        if isinstance(files, str):
            self._files = [files]
        elif isinstance(files, Iterable):  # test items of iterable
            files = list(files)
            if all(isinstance(f, str) for f in files):
                self._files = files
            else:
                raise TypeError("All passed files must be type str")
        else:
            raise TypeError("Files must be a string filename or a list of such names.")

        # set open mode:
        if mode not in {"r", "rb"}:
            raise ValueError("Mode must be one of 'r', 'rb'")
        self._mode = mode

        if isinstance(header_comment_char, str) and mode == "rb":
            self._header_comment_char = header_comment_char.encode()
        else:
            self._header_comment_char = header_comment_char

    @property
    def filenames(self) -> List[str]:
        return self._files

    def __len__(self):
        """Return the length of the Reader object.

        Notes
        -----
        This function requires reading the complete file, and should typically not be
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
    def size(self) -> int:
        """return the collective size of all files being read in bytes"""
        return sum(os.stat(f).st_size for f in self._files)

    def select_record_indices(self, indices: Set) -> Generator:
        """Iterate over provided indices only, skipping other records.

        Parameters
        ----------
        indices : Set[int]
            indices to include in the output

        Yields
        ------
        record, str
            records from file corresponding to indices

        """
        indices = copy(
            indices
        )  # passed indices is a reference, need own copy to modify
        for idx, record in enumerate(self):
            if idx in indices:
                yield record
                indices.remove(idx)

                # stopping condition
                if not indices:
                    break


def zip_readers(*readers, indices=None) -> Generator:
    """Zip together multiple reader objects, yielding records simultaneously.

    If indices is passed, only return lines in file that correspond to indices

    Parameters
    ----------
    *readers : List[Reader]
        Reader objects to simultaneously iterate over
    indices : Set[int], optional
        indices to include in the output

    Yields
    ------
    records : Tuple[str]
        one record per reader passed

    """
    if indices:
        iterators = zip(*(r.select_record_indices(indices) for r in readers))
    else:
        iterators = zip(*readers)
    for record_tuple in iterators:
        yield record_tuple
