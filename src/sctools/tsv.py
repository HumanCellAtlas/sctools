"""
Tools for Manipulating TSV format files
===========================================

.. currentmodule:: sctools

This module provides functions and classes to subsample lines from tsv files that correspond to
alignments, split tsv files into chunks, and iterate over sorted tsv files by one or more tags


Methods
-------
iter_tag_groups_from_tsv                         function to iterate over lines by an arbitrary tag pos


"""

from typing import Iterator, Generator
from io import TextIOWrapper

def iter_tag_groups_from_tsv(
    tsv_iterator: Iterator[TextIOWrapper], filter_null: bool = False
) -> Generator:
    """Iterates over lines and yields them grouped by the provided tag value

    Parameters
    ----------
    tag : str
        BAM tag to group over
    tsv_iterator : Iterator[pysam.AlignedSegment]
        open tsv file that can be iterated over
    filter_null : bool, optional
        If False, all lines that lack the requested tag are yielded together. Else, all lines
        that lack the tag will be discarded (default = False).

    Yields
    ------
    grouped_by_tag : Iterator[pysam.AlignedSegment]
        lines sharing a unique value of tag
    current_tag : str
        the tag that lines in the group all share

    """

    # get first read and tag set
    lines = [next(tsv_iterator)]
    fields = [x.strip() for x in lines[0].split("\t")]
    current_tag = fields[0]

    # now iterate over alignment sets
    for alignment in tsv_iterator:
        fields = [x.strip() for x in alignment.split("\t")]
        next_tag = fields[0]

        if next_tag == current_tag:
            lines.append(alignment)
        else:
            # only yield if the tag is non-null or filter_null is false
            if not filter_null or current_tag is not None:
                yield lines, current_tag
            # reset to next group
            lines = [alignment]
            current_tag = next_tag

    if not filter_null or current_tag is not None:
        yield lines, current_tag
