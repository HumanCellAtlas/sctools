"""
Tools for Manipulating SAM/BAM format files
===========================================

.. currentmodule:: sctools

This module provides functions and classes to subsample reads from bam files that correspond to
specific chromosomes, split bam files into chunks, assign tags to bam files from paired fastq
records, and iterate over sorted bam files by one or more tags

This module makes heavy use of the pysam wrapper for HTSlib, a high-performance c-library designed
to manipulate sam files

Methods
-------
iter_tag_groups                         function to iterate over reads by an arbitrary tag
iter_cell_barcodes                      wrapper for iter_tag_groups that iterates over cell barcode tags
iter_genes                              wrapper for iter_tag_groups that iterates over gene tags
iter_molecules                          wrapper for iter_tag_groups that iterates over molecule tags
sort_by_tags_and_queryname              sort bam by given list of zero or more tags, followed by query name
verify_sort                             verifies whether bam is correctly sorted by given list of tags, then query name

Classes
-------
SubsetAlignments                        class to extract reads specific to requested chromosome(s)
Tagger                                  class to add tags to sam/bam records from paired fastq records
AlignmentSortOrder                      abstract class to represent alignment sort orders
QueryNameSortOrder                      alignment sort order by query name
TagSortableRecord                       class to facilitate sorting of pysam.AlignedSegments
SortError                               error raised when sorting is incorrect

References
----------
htslib : https://github.com/samtools/htslib

"""

import functools
from functools import partial
from functools import reduce
import math
import os
import warnings
from abc import abstractmethod
from typing import (
    Iterator,
    Iterable,
    Generator,
    List,
    Set,
    Dict,
    Union,
    Tuple,
    Callable,
    Any,
    Optional,
)

import pysam
import shutil
import multiprocessing
import random

from . import consts

# File descriptor to write log messages to
STDERR = 2


class SubsetAlignments:
    """Wrapper for pysam/htslib that extracts reads corresponding to requested chromosome(s)

    Parameters
    ----------
    alignment_file : str
        sam or bam file
    open_mode : {'r', 'rb', None}, optional
        open mode for pysam.AlignmentFile. 'r' indicates a sam file, 'rb' indicates a bam file,
        and None attempts to autodetect based on the file suffix (Default = None)

    Methods
    -------
    indices_by_chromosome
        returns indices to line numbers containing the requested number of reads for a specified
        chromosome

    Notes
    -----
    samtools is a good general-purpose tool for that is capable of most subsampling tasks. It is a
    good idea to check the samtools documentation when approaching these types of tasks.

    References
    ----------
    samtools documentation : http://www.htslib.org/doc/samtools.html

    """

    def __init__(self, alignment_file: str, open_mode: str = None):
        if open_mode is None:
            if alignment_file.endswith('.bam'):
                open_mode = 'rb'
            elif alignment_file.endswith('.sam'):
                open_mode = 'r'
            else:
                raise ValueError(
                    f'Could not autodetect file type for alignment_file {alignment_file} (detectable suffixes: '
                    f'.sam, .bam)'
                )
        self._file: str = alignment_file
        self._open_mode: str = open_mode

    def indices_by_chromosome(
        self, n_specific: int, chromosome: str, include_other: int = 0
    ) -> Union[List[int], Tuple[List[int], List[int]]]:
        """Return the list of first `n_specific` indices of reads aligned to `chromosome`.

        Parameters
        ----------
        n_specific : int
            Number of aligned reads to return indices for
        chromosome : str
            Only reads from this chromosome are considered valid
        include_other : int, optional
            The number of reads to include that are NOT aligned to chromosome. These can be aligned
            or unaligned reads (default = 0).

        Returns
        -------
        chromosome_indices : List[int]
            list of indices to reads aligning to `chromosome`
        other_indices : List[int], optional
            list of indices to reads NOT aligning to chromosome, only returned if include_other is
            not 0.

        """

        # acceptable chromosomes
        valid_chromosomes = [str(i) for i in range(1, 23)] + ['M', 'MT', 'X', 'Y']
        valid_chromosomes.extend(['chr' + v for v in valid_chromosomes])

        # check chromosome
        if isinstance(chromosome, int) and chromosome < 23:
            chromosome = str(chromosome)  # try to convert
        if chromosome not in valid_chromosomes:
            warnings.warn(
                'chromsome %s not in list of expected chromosomes: %r'
                % (chromosome, valid_chromosomes)
            )

        with pysam.AlignmentFile(self._file, self._open_mode) as fin:
            chromosome = str(chromosome)
            chromosome_indices = []
            other_indices = []

            for i, record in enumerate(fin):

                if not record.is_unmapped:  # record is mapped
                    if chromosome == record.reference_name:
                        if len(chromosome_indices) < n_specific:
                            chromosome_indices.append(i)
                    elif len(other_indices) < include_other:
                        other_indices.append(i)
                elif len(other_indices) < include_other:  # record is not mapped
                    other_indices.append(i)

                # check termination condition (we have the requisite number of reads
                if (
                    len(chromosome_indices) == n_specific
                    and len(other_indices) == include_other
                ):
                    break

        if len(chromosome_indices) < n_specific or len(other_indices) < include_other:
            warnings.warn(
                'Only %d unaligned and %d reads aligned to chromosome %s were found in'
                '%s'
                % (len(other_indices), len(chromosome_indices), chromosome, self._file)
            )

        if include_other != 0:
            return chromosome_indices, other_indices
        else:
            return chromosome_indices


class Tagger:
    """Add tags to a bam file from tag generators.

    Parameters
    ----------
    bam_file : str
        Bam file that tags are to be added to.

    Methods
    -------
    tag
        tag bam records given tag_generators (often generated from paired bam or fastq files)
        # todo this should probably be wrapped up in __init__ to make this more function-like
    """

    def __init__(self, bam_file: str) -> None:
        if not isinstance(bam_file, str):
            raise TypeError(
                f'The argument "bam_file" must be of type str, not {type(bam_file)}'
            )
        self.bam_file = bam_file

    # todo add type to tag_generators (make sure it doesn't introduce import issues
    def tag(self, output_bam_name: str, tag_generators) -> None:
        """Add tags to bam_file.

        Given a bam file and tag generators derived from files sharing the same sort order,
        adds tags to the .bam file, and writes the resulting file to output_bam_name.

        Parameters
        ----------
        output_bam_name : str
            Name of output tagged bam.
        tag_generators : List[fastq.TagGenerator]
            list of generators that yield fastq.Tag objects

        """
        with pysam.AlignmentFile(
            self.bam_file, 'rb', check_sq=False
        ) as inbam, pysam.AlignmentFile(
            output_bam_name, 'wb', template=inbam
        ) as outbam:

            # zip up all the iterators
            for *tag_sets, sam_record in zip(*tag_generators, inbam):
                for tag_set in tag_sets:
                    for tag in tag_set:
                        sam_record.set_tag(*tag)
                outbam.write(sam_record)


def get_barcodes_from_bam(
    in_bam: str, tags: List[str], raise_missing: bool
) -> Set[str]:
    """ Get all the distinct barcodes from a bam

    :param in_bam: str
        Input bam file.
    :param tags: List[str]
        Tags in the bam that might contain barcodes.
    :param raise_missing: bool
        Raise an error if no barcodes can be found.
    :return: set
        A set of barcodes found in the bam
        This set will not contain a None value
    """
    barcodes = set()
    # Get all the Barcodes from the BAM
    with pysam.AlignmentFile(in_bam, 'rb', check_sq=False) as input_alignments:
        for alignment in input_alignments:
            barcode = get_barcode_for_alignment(alignment, tags, raise_missing)
            # If no provided tag was found on the record that had a non-null value
            if barcode is not None:
                barcodes.add(barcode)
    return barcodes


def get_barcode_for_alignment(
    alignment: pysam.AlignedSegment, tags: List[str], raise_missing: bool
) -> str:
    """ Get the barcode for an Alignment

    :param alignment: pysam.AlignedSegment
        An Alignment from pysam.
    :param tags: List[str]
        Tags in the bam that might contain barcodes.
    :param raise_missing: bool
        Raise an error if no barcodes can be found.
    :return: str
        A barcode for the alignment, or None if one is not found and raise_missing is False.
    """
    alignment_barcode = None
    for tag in tags:
        # The non-existent barcode should be the exceptional case, so try/except is faster than if/else
        try:
            alignment_barcode = alignment.get_tag(tag)
            break  # Got the key, don't bother getting the next tag
        except KeyError:
            continue  # Try to get the next tag

    if raise_missing and alignment_barcode is None:
        raise RuntimeError(
            'Alignment encountered that is missing {} tag(s).'.format(tags)
        )

    return alignment_barcode


def write_barcodes_to_bins(
    in_bam: str, tags: List[str], barcodes_to_bins: Dict[str, int], raise_missing: bool
) -> List[str]:
    """ Write barcodes to appropriate shards as defined by barcodes_to_bins

    :param in_bam: str
        The bam file to read.
    :param tags: List[str]
        Tags in the bam that might contain barcodes.
    :param barcodes_to_bins: Dict[str, int]
        A Dict from barcode to bin. All barcodes of the same type need to be written to the same shard.
        These numbered shards are merged after parallelization so that all alignments with the same
        barcode are in the same bam.
    :param raise_missing: bool
        Raise an error if no barcodes can be found.
    :return: A list of paths to the written shards.
    """
    # Create all the output files
    with pysam.AlignmentFile(in_bam, 'rb', check_sq=False) as input_alignments:

        # We need a random int appended to the dirname to make sure input bams with the same name don't clash
        dirname = (
            os.path.splitext(os.path.basename(in_bam))[0]
            + "_"
            + str(random.randint(0, 10000))
        )
        os.makedirs(dirname)

        files = []
        bins = list(set(barcodes_to_bins.values()))
        filepaths = []
        for i in range(len(bins)):
            out_bam_name = os.path.realpath(dirname) + ("/" + dirname + '_%d.bam' % i)
            filepaths.append(out_bam_name)
            # For now, bam writing uses one thread for compression. Better logic could support more threads without
            # starving the machine for resources
            open_bam = pysam.AlignmentFile(
                out_bam_name, 'wb', template=input_alignments
            )
            files.append(open_bam)

        # Loop over input; check each tag in priority order and partition barcodes into files based
        # on the highest priority tag that is identified
        for alignment in input_alignments:
            barcode = get_barcode_for_alignment(alignment, tags, raise_missing)
            if barcode is not None:
                # Find or set the file associated with the tag and write the record to the correct file
                out_file = files[barcodes_to_bins[barcode]]
                out_file.write(alignment)

    for file in files:
        file.close()

    return filepaths


def merge_bams(bams: List[str]) -> str:
    """ Merge input bams using samtools.

    This cannot be a local function within `split` because then Python "cannot pickle a local object".
    :param bams: Bams to merge.
    :return: The output bam name.
    """
    bam_name = os.path.realpath(bams[0] + ".bam")
    bams_to_merge = bams[1:]
    pysam.merge('-c', '-p', bam_name, *bams_to_merge)
    return bam_name


def split(
    in_bams: List[str],
    out_prefix: str,
    tags: List[str],
    approx_mb_per_split: float = 1000,
    raise_missing: bool = True,
    num_threads: int = None,
) -> List[str]:
    """split `in_bam` by tag into files of `approx_mb_per_split`

    Parameters
    ----------
    in_bams : str
        Input bam files.
    out_prefix : str
        Prefix for all output files; output will be named as prefix_n where n is an integer equal
        to the chunk number.
    tags : List[str]
        The bam tags to split on. The tags are checked in order, and sorting is done based on the
        first identified tag. Further tags are only checked if the first tag is missing. This is
        useful in cases where sorting is executed over a corrected barcode, but some records only
        have a raw barcode.
    approx_mb_per_split : float
        The target file size for each chunk in mb
    raise_missing : bool, optional
        if True, raise a RuntimeError if a record is encountered without a tag. Else silently
        discard the record (default = True)
    num_threads : int, optional
        The number of threads to parallelize over. If not set, will use all available threads.

    Returns
    -------
    output_filenames : List[str]
        list of filenames of bam chunks

    Raises
    ------
    ValueError
        when `tags` is empty
    RuntimeError
        when `raise_missing` is true and any passed read contains no `tags`

    """

    if len(tags) == 0:
        raise ValueError('At least one tag must be passed')

    if num_threads is None:
        num_threads = multiprocessing.cpu_count()

    # find correct number of subfiles to spawn
    bam_mb = sum(map(lambda in_bam: os.path.getsize(in_bam) * 1e-6, in_bams))
    n_subfiles = int(math.ceil(bam_mb / approx_mb_per_split))
    if n_subfiles > consts.MAX_BAM_SPLIT_SUBFILES_TO_WARN:
        warnings.warn(
            f'Number of requested subfiles ({n_subfiles}) exceeds '
            f'{consts.MAX_BAM_SPLIT_SUBFILES_TO_WARN}; this may cause OS errors by exceeding fid limits'
        )
    if n_subfiles > consts.MAX_BAM_SPLIT_SUBFILES_TO_RAISE:
        raise ValueError(
            f'Number of requested subfiles ({n_subfiles}) exceeds '
            f'{consts.MAX_BAM_SPLIT_SUBFILES_TO_RAISE}; this will usually cause OS errors, '
            f'think about increasing max_mb_per_split.'
        )

    full_pool = multiprocessing.Pool(num_threads)

    # Get all the barcodes over all the bams
    os.write(STDERR, b'Retrieving barcodes from bams\n')
    result = full_pool.map(
        partial(get_barcodes_from_bam, tags=tags, raise_missing=raise_missing), in_bams
    )

    barcodes_list = list(reduce(lambda set1, set2: set1.union(set2), result))
    os.write(STDERR, b'Retrieved barcodes from bams\n')

    # Create the barcodes to bin mapping
    os.write(STDERR, b'Allocating bins\n')
    barcodes_to_bins_dict = {}

    # barcodes_list will always contain non-None elements from get_barcodes_from_bam
    if len(barcodes_list) <= n_subfiles:
        for barcode_index in range(len(barcodes_list)):
            barcodes_to_bins_dict[barcodes_list[barcode_index]] = barcode_index
    else:
        for barcode_index in range(len(barcodes_list)):
            file_index = barcode_index % n_subfiles
            barcodes_to_bins_dict[barcodes_list[barcode_index]] = file_index

    # Split the bams by barcode in parallel
    os.write(STDERR, b'Splitting the bams by barcode\n')
    # Samtools needs a thread for compression, so we leave half the given threads open.
    write_pool_threads = math.ceil(num_threads / 2) if num_threads > 2 else 1
    write_pool = multiprocessing.Pool(write_pool_threads)
    scattered_split_result = write_pool.map(
        partial(
            write_barcodes_to_bins,
            tags=list(tags),
            raise_missing=raise_missing,
            barcodes_to_bins=barcodes_to_bins_dict,
        ),
        in_bams,
    )

    bin_indices = list(set(barcodes_to_bins_dict.values()))
    bins = list(map(lambda index: ["{}_{}".format(out_prefix, index)], bin_indices))

    for shard_index in range(len(scattered_split_result)):
        shard = scattered_split_result[shard_index]
        for file_index in range(len(shard)):
            bins[file_index].append(shard[file_index])

    write_pool.close()

    # Recombine the binned bams
    os.write(STDERR, b'Merging temporary bam files\n')
    merged_bams = full_pool.map(partial(merge_bams), bins)

    os.write(STDERR, b'deleting temporary files\n')
    for paths in scattered_split_result:
        shutil.rmtree(os.path.dirname(paths[0]))

    full_pool.close()

    return merged_bams


# todo change this to throw away "None" reads instead of appending them if we are filtering them
def iter_tag_groups(
    tag: str, bam_iterator: Iterator[pysam.AlignedSegment], filter_null: bool = False
) -> Generator:
    """Iterates over reads and yields them grouped by the provided tag value

    Parameters
    ----------
    tag : str
        BAM tag to group over
    bam_iterator : Iterator[pysam.AlignedSegment]
        open bam file that can be iterated over
    filter_null : bool, optional
        If False, all reads that lack the requested tag are yielded together. Else, all reads
        that lack the tag will be discarded (default = False).

    Yields
    ------
    grouped_by_tag : Iterator[pysam.AlignedSegment]
        reads sharing a unique value of tag
    current_tag : str
        the tag that reads in the group all share

    """

    # get first read and tag set
    reads = [next(bam_iterator)]
    try:
        current_tag = reads[0].get_tag(tag)
    except KeyError:
        current_tag = None  # null tag is a category that gets emitted

    # now iterate over alignment sets
    for alignment in bam_iterator:
        try:
            next_tag = alignment.get_tag(tag)
        except KeyError:
            next_tag = None  # null tag is a category that we will emit
        if next_tag == current_tag:
            reads.append(alignment)
        else:
            # only yield if the tag is non-null or filter_null is false
            if not filter_null or current_tag is not None:
                yield iter(reads), current_tag
            # reset to next group
            reads = [alignment]
            current_tag = next_tag

    if not filter_null or current_tag is not None:
        yield iter(reads), current_tag


def iter_molecule_barcodes(bam_iterator: Iterator[pysam.AlignedSegment]) -> Generator:
    """Iterate over all the molecules of a bam file sorted by molecule.

    Parameters
    ----------
    bam_iterator : Iterator[pysam.AlignedSegment]
        open bam file that can be iterated over

    Yields
    ------
    grouped_by_tag : Iterator[pysam.AlignedSegment]
        reads sharing a unique molecule barcode tag
    current_tag : str
        the molecule barcode that records in the group all share

    """
    return iter_tag_groups(
        tag=consts.MOLECULE_BARCODE_TAG_KEY, bam_iterator=bam_iterator
    )


def iter_cell_barcodes(bam_iterator: Iterator[pysam.AlignedSegment]) -> Generator:
    """Iterate over all the cells of a bam file sorted by cell.

    Parameters
    ----------
    bam_iterator : Iterator[pysam.AlignedSegment]
        open bam file that can be iterated over

    Yields
    ------
    grouped_by_tag : Iterator[pysam.AlignedSegment]
        reads sharing a unique cell barcode tag
    current_tag : str
        the cell barcode that reads in the group all share

    """
    return iter_tag_groups(tag=consts.CELL_BARCODE_TAG_KEY, bam_iterator=bam_iterator)


def iter_genes(bam_iterator: Iterator[pysam.AlignedSegment]) -> Generator:
    """Iterate over all the cells of a bam file sorted by gene.

    Parameters
    ----------
    bam_iterator : Iterator[pysam.AlignedSegment]
        open bam file that can be iterated over

    Yields
    ------
    grouped_by_tag : Iterator[pysam.AlignedSegment]
        reads sharing a unique gene name tag
    current_tag : str
        the gene id that reads in the group all share

    """
    return iter_tag_groups(tag=consts.GENE_NAME_TAG_KEY, bam_iterator=bam_iterator)


def get_tag_or_default(
    alignment: pysam.AlignedSegment, tag_key: str, default: Optional[str] = None
) -> Optional[str]:
    """Extracts the value associated to `tag_key` from `alignment`, and returns a default value
    if the tag is not present."""
    try:
        return alignment.get_tag(tag_key)
    except KeyError:
        return default


class AlignmentSortOrder:
    """The base class of alignment sort orders."""

    @property
    @abstractmethod
    def key_generator(self) -> Callable[[pysam.AlignedSegment], Any]:
        """Returns a callable function that calculates a sort key from given pysam.AlignedSegment."""
        raise NotImplementedError


class QueryNameSortOrder(AlignmentSortOrder):
    """Alignment record sort order by query name."""

    @staticmethod
    def get_sort_key(alignment: pysam.AlignedSegment) -> str:
        return alignment.query_name

    @property
    def key_generator(self):
        return QueryNameSortOrder.get_sort_key

    def __repr__(self) -> str:
        return 'query_name'


@functools.total_ordering
class TagSortableRecord(object):
    """Wrapper for pysam.AlignedSegment that facilitates sorting by tags and query name."""

    def __init__(
        self,
        tag_keys: Iterable[str],
        tag_values: Iterable[str],
        query_name: str,
        record: pysam.AlignedSegment = None,
    ) -> None:
        self.tag_keys = tag_keys
        self.tag_values = tag_values
        self.query_name = query_name
        self.record = record

    @classmethod
    def from_aligned_segment(
        cls, record: pysam.AlignedSegment, tag_keys: Iterable[str]
    ) -> 'TagSortableRecord':
        """Create a TagSortableRecord from a pysam.AlignedSegment and list of tag keys"""
        assert record is not None
        tag_values = [get_tag_or_default(record, key, '') for key in tag_keys]
        query_name = record.query_name
        return cls(tag_keys, tag_values, query_name, record)

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, TagSortableRecord):
            return NotImplemented
        self.__verify_tag_keys_match(other)
        for (self_tag_value, other_tag_value) in zip(self.tag_values, other.tag_values):
            if self_tag_value < other_tag_value:
                return True
            elif self_tag_value > other_tag_value:
                return False
        return self.query_name < other.query_name

    def __eq__(self, other: object) -> bool:
        # TODO: Add more error checking
        if not isinstance(other, TagSortableRecord):
            return NotImplemented
        self.__verify_tag_keys_match(other)
        for (self_tag_value, other_tag_value) in zip(self.tag_values, other.tag_values):
            if self_tag_value != other_tag_value:
                return False
        return self.query_name == other.query_name

    def __verify_tag_keys_match(self, other) -> None:
        if self.tag_keys != other.tag_keys:
            format_str = 'Cannot compare records using different tag lists: {0}, {1}'
            raise ValueError(format_str.format(self.tag_keys, other.tag_keys))

    def __str__(self) -> str:
        return self.__repr__()

    def __repr__(self) -> str:
        format_str = 'TagSortableRecord(tags: {0}, tag_values: {1}, query_name: {2}'
        return format_str.format(self.tag_keys, self.tag_values, self.query_name)


def sort_by_tags_and_queryname(
    records: Iterable[pysam.AlignedSegment], tag_keys: Iterable[str]
) -> Iterable[pysam.AlignedSegment]:
    """Sorts the given bam records by the given tags, followed by query name.
    If no tags are given, just sorts by query name.
    """
    tag_sortable_records = (
        TagSortableRecord.from_aligned_segment(r, tag_keys) for r in records
    )
    sorted_records = sorted(tag_sortable_records)
    aligned_segments = (r.record for r in sorted_records)
    return aligned_segments


def verify_sort(records: Iterable[TagSortableRecord], tag_keys: Iterable[str]) -> None:
    """Raise AssertionError if the given records are not correctly sorted by the given tags and query name"""
    # Setting tag values and query name to empty string ensures first record will never be less than old_record
    old_record = TagSortableRecord(
        tag_keys=tag_keys, tag_values=['' for _ in tag_keys], query_name='', record=None
    )
    i = 0
    for record in records:
        i += 1
        if not record >= old_record:
            msg = 'Records {0} and {1} are not in correct order:\n{1}:{2} \nis less than \n{0}:{3}'
            raise SortError(msg.format(i - 1, i, record, old_record))
        old_record = record


class SortError(Exception):
    pass
