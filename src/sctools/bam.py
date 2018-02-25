import warnings
import os
import math
import pysam
from itertools import cycle
from typing import Iterator, Generator

"""
unlike fastq and gtf which lack flexible iterators, the pysam wrapper provides an excellent iterator 
over sam records that would be difficult to improve upon in pure python. 

This module provides a few additional methods that extend the functionality of pysam
"""


class SubsetAlignments:
    """Wrapper for pysam/htslib that enables non-standard filtering of alignments in a bamfile

    :method indices_by_chromosome: extract reads specific to a chromosome
    """

    def __init__(self, alignment_file, open_mode=None):
        """
        :param str alignment_file: sam or bam file.
        :param str open_mode: optional, mode to read file. Will be autodetected by file type if
          the file contains the correct suffix for its type.
        """

        if open_mode is None:
            if alignment_file.endswith('.bam'):
                open_mode = 'rb'
            elif alignment_file.endswith('.sam'):
                open_mode = 'r'
            else:
                raise ValueError('could not autodetect file type for alignment_file %s '
                                 '(detectible suffixes: .sam, .bam)' % alignment_file)
        self._file = alignment_file
        self._open_mode = open_mode

    def indices_by_chromosome(self, n_specific, chromosome, include_other=0):
        """Return the list of first n_specific indices of reads aligned to selected chromosome.

        If desired, will also return non-specific indices in a second list (can serve as negative
        control reads).

        :param int n_specific: number of aligned reads to return indices for
        :param str chromosome: only reads from this chromosome are considered valid
        :param int include_other: optional, (default=0), the number of reads to include that are
          NOT aligned to chromosome (could be aligned or unaligned read)
        :return [int]: list of indices to reads aligning to chromosome
        :return [int]: list of indices to reads NOT aligning to chromosome, only returned if
          include_other is not 0.
        """

        # acceptable chromosomes
        valid_chromosomes = [str(i) for i in range(1, 23)] + ['M', 'MT', 'X', 'Y']
        valid_chromosomes.extend(['chr' + v for v in valid_chromosomes])

        # check chromosome
        if isinstance(chromosome, int) and chromosome < 23:
            chromosome = str(chromosome)  # try to convert
        if chromosome not in valid_chromosomes:
            warnings.warn('chromsome %s not in list of expected chromosomes: %r' %
                          (chromosome, valid_chromosomes))

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
                if len(chromosome_indices) == n_specific and len(other_indices) == include_other:
                    break

        if len(chromosome_indices) < n_specific or len(other_indices) < include_other:
            warnings.warn('Only %d unaligned and %d reads aligned to chromosome %s were found in' 
                          '%s' % (len(other_indices), len(chromosome_indices),
                                  chromosome, self._file))

        if include_other != 0:
            return chromosome_indices, other_indices
        else:
            return chromosome_indices


class Tagger:
    """Add tags to a bam file from tag generators

    :method tag: add tags from an arbitrary number of tag generators. Typical use is to extract
      barcodes from fastq files.
    """

    def __init__(self, bam_file):
        """
        :param str bam_file: location of bam file
        """
        if not isinstance(bam_file, str):
            raise TypeError('bam_file must be type str, not %s' % type(bam_file))
        self.bam_file = bam_file

    def tag(self, output_bam_name, tag_generators):
        """
        given a bam file and tag generators derived from files sharing the same sort order,
        adds tags to the .bam file, writes the resulting file to output_bam_name.

        :param str output_bam_name: name of output tagged bam.
        :param [fastq.TagGenerator] tag_generators: generators that yield Tag objects
          (see fastq.Tag)
        """
        with pysam.AlignmentFile(self.bam_file, 'rb', check_sq=False) as inbam, \
                pysam.AlignmentFile(output_bam_name, 'wb', template=inbam) as outbam:

            # zip up all the iterators
            for *tag_sets, sam_record in zip(*tag_generators, inbam):
                for tag_set in tag_sets:
                    for tag in tag_set:
                        sam_record.set_tag(*tag)
                outbam.write(sam_record)


def split(in_bam, out_prefix, tag, approx_mb_per_split=1000, raise_missing=True):
    """
    split a bam file into subfiles by tag, calculating the number of splits so that the chunks are
    approximately approx_mb_per_split

    :param str in_bam: input bam file
    :param str out_prefix:  prefix for all output files; output will be named as prefix_n where n
      is an integer equal to the chunk number.
    :param tag: the bam tag to split on
    :param float approx_mb_per_split: the target file size for each chunk in mb
    :param raise_missing: default=True, if True, raise a RuntimeError if a record is encountered
      without a tag. Else silently discard the record.


    :return [str]: output file names
    """

    def _cleanup(files_to_counts, files_to_names, rm_all=False):
        """close files, remove any empty files.

        :param dict files_to_counts:
        :param dict files_to_names:
        :param bool rm_all: indicates all files should be removed, regardless of count number.
        :return:
        """
        for bamfile, count in files_to_counts.items():
            # corner case: clean up files that were created, but didn't get data because
            # n_cell < n_barcode
            bamfile.close()
            if count == 0 or rm_all:
                os.remove(files_to_names[bamfile])
                del files_to_names[bamfile]

    # find correct number of subfiles to spawn
    bam_mb = os.path.getsize(in_bam) * 1e-6
    n_subfiles = int(math.ceil(bam_mb / approx_mb_per_split))
    if n_subfiles > 500:
        warnings.warn('Number of requested subfiles (%d) exceeds 500; this may cause OS errors by '
                      'exceeding fid limits' % n_subfiles)
    if n_subfiles > 1000:
        raise ValueError('Number of requested subfiles (%d) exceeds 1000; this will usually cause '
                         'OS errors, please increase max_mb_per_split.' % n_subfiles)

    # create all the output files
    with pysam.AlignmentFile(in_bam, 'rb', check_sq=False) as input_alignments:

        # map files to counts
        files_to_counts = {}
        files_to_names = {}
        for i in range(n_subfiles):
            out_bam_name = out_prefix + '_%d.bam' % i
            open_bam = pysam.AlignmentFile(out_bam_name, 'wb', template=input_alignments)
            files_to_counts[open_bam] = 0
            files_to_names[open_bam] = out_bam_name

        # cycler over files to assign new barcodes to next file
        file_cycler = cycle(files_to_counts.keys())

        # create an empty map for barcodes to files
        tags_to_files = {}

        # loop over input, partitioning barcodes into files
        for alignment in input_alignments:
            try:
                tag_content = alignment.get_tag(tag)
            except KeyError:
                if raise_missing:
                    _cleanup(files_to_counts, files_to_names, rm_all=True)
                    raise RuntimeError(
                        'alignment encountered that is missing %s tag.' % tag)
                else:
                    continue
            try:
                out_file = tags_to_files[tag_content]
            except KeyError:
                out_file = next(file_cycler)
                tags_to_files[tag_content] = out_file
            out_file.write(alignment)
            files_to_counts[out_file] += 1

    _cleanup(files_to_counts, files_to_names)
    return list(files_to_names.values())

def iter_tag_groups(tag: str, bam_iterator: Iterator[pysam.AlignedSegment]) -> Generator:
    # get first read and tag set
    reads = []
    while not reads:

        try:
            read = next(bam_iterator)
            current_tag = read.get_tag(tag)
            reads.append(read)
            break
        except KeyError:
            continue
        except StopIteration:
            raise ValueError('No reads contained tag %s' % tag)

    for alignment in bam_iterator:
        try:
            next_tag = alignment.get_tag(tag)
        except KeyError:
            continue
        if next_tag == current_tag:
            reads.append(alignment)
        else:
            yield iter(reads), current_tag
            reads = [alignment]
            current_tag = next_tag
    yield iter(reads), current_tag


def iter_molecule_barcodes(bam_iterator: Iterator[pysam.AlignedSegment]) -> Generator:
    """
    function to iterate over all the molecules of a bam file sorted by molecule
    """
    return iter_tag_groups(tag='UB', bam_iterator=bam_iterator)


def iter_cell_barcodes(bam_iterator: Iterator[pysam.AlignedSegment]) -> Generator:
    """
    function to iterate over all the cells of a bam file sorted by cell
    """
    return iter_tag_groups(tag='CB', bam_iterator=bam_iterator)


def iter_genes(bam_iterator: Iterator[pysam.AlignedSegment]) -> Generator:
    """
    function to iterate over all the cells of a bam file sorted by gene
    """
    return iter_tag_groups(tag='GE', bam_iterator=bam_iterator)
