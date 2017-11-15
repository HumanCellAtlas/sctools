import warnings
import pysam

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
                pysam.AlignmentFile(output_bam_name, 'wb', header=inbam.header) as outbam:

            # zip up all the iterators
            for *tag_sets, sam_record in zip(*tag_generators, inbam):
                for tag_set in tag_sets:
                    for tag in tag_set:
                        sam_record.set_tag(*tag)
                outbam.write(sam_record)
