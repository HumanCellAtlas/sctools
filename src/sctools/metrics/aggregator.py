"""
Sequence Metric Aggregators
===========================

.. currentmodule:: sctools.metrics

This module provides classes useful for aggregating metric information for individual cells or
genes. These classes consume BAM files that have been pre-sorted such that all sequencing reads
that correspond to the molecules of a cell (CellMetrics) or the molecules of a gene (GeneMetrics)
are yielded sequentially.

Classes
-------

.. autosummary::
   :toctree: generated/

   MetricAggregatorBase         Aggregator Base Class
   GeneMetrics                  Class to iteratively calculate metrics for a gene (by molecule)
   CellMetrics                  Class to iteratively calculate metrics for a cell (by molecule)

Notes
-----
This module can be rewritten with dataclass when python 3.7 stabilizes, see
https://www.python.org/dev/peps/pep-0557/


See Also
--------
sctools.metrics.gatherer
sctools.metrics.merge
sctools.metrics.writer

"""

from collections import Counter
from typing import Iterable, Tuple, Counter, List, Sequence

import numpy as np
import pysam

from sctools import consts
from sctools.stats import OnlineGaussianSufficientStatistic


class MetricAggregator:
    """Metric Aggregator Base Class

    The ``MetricAggregator`` class defines a set of metrics that can be extracted from an
    aligned bam file. It defines all the metrics that are general across genes and cells. This
    class is subclassed by ``GeneMetrics`` and ``CellMetrics``, which define data-specific metrics
    in the ``parse_extra_fields`` method. An instance of ``GeneMetrics`` or ``CellMetrics`` is
    instantiated for each gene or molecule in a bam file, respectively.

    Attributes
    ----------
    n_reads : int
        The number of reads associated with this entity
    noise_reads : int, NotImplemented
        Number of reads that are categorized by 10x genomics cellranger as "noise". Refers to
        long polymers, or reads with high numbers of N (ambiguous) nucleotides
    perfect_molecule_barcodes : int
        The number of reads with molecule barcodes that have no errors (cell barcode tag == raw barcode tag)
    reads_mapped_exonic : int
        The number of reads for this entity that are mapped to exons
    reads_mapped_intronic : int
        The number of reads for this entity that are mapped to introns
    reads_mapped_utr : int
        The number of reads for this entity that are mapped to 3' untranslated regions (UTRs)
    reads_mapped_uniquely : int
        The number of reads mapped to a single unambiguous location in the genome
    reads_mapped_multiple : int
        The number of reads mapped to multiple genomic positions with equal confidence
        # todo make sure equal confidence is accurate
    duplicate_reads : int
        The number of reads that are duplicates (see README.md for defition of a duplicate)
    spliced_reads : int
        The number of reads that overlap splicing junctions
    antisense_reads : int
        The number of reads that are mapped to the antisense strand instead of the transcribed
        strand
    molecule_barcode_fraction_bases_above_30_mean : float
        The average fraction of bases in molecule barcodes that receive quality scores greater than
        30 across the reads of this entity
    molecule_barcode_fraction_bases_above_30_variance : float
        The variance in the fraction of bases in molecule barcodes that receive quality scores
        greater than 30 across the reads of this entity
    genomic_reads_fraction_bases_quality_above_30_mean : float
        The average fraction of bases in the genomic read that receive quality scores greater than
        30 across the reads of this entity (included for 10x cell ranger count comparison)
    genomic_reads_fraction_bases_quality_above_30_variance : float
        The variance in the fraction of bases in the genomic read that receive quality scores
        greater than 30 across the reads of this entity (included for 10x cell ranger count
        comparison)
    genomic_read_quality_mean : float
        Average quality of Illumina base calls in the genomic reads corresponding to this entity
    genomic_read_quality_variance : float
        Variance in quality of Illumina base calls in the genomic reads corresponding to this
        entity
    n_molecules : float
        Number of molecules corresponding to this entity. See README.md for the definition of a
        Molecule
    n_fragments : float
        Number of fragments corresponding to this entity. See README.md for the definition of a
        Fragment
    reads_per_molecule : float
        The average number of reads associated with each molecule in this entity
    reads_per_fragment : float
        The average number of reads associated with each fragment in this entity
    fragments_per_molecule : float
        The average number of fragments associated with each molecule in this entity
    fragments_with_single_read_evidence : int
        The number of fragments associated with this entity that are observed by only one read
    molecules_with_single_read_evidence : int
        The number of molecules associated with this entity that are observed by only one read

    Methods
    -------
    parse_extra_fields(tags, record), NotImplemented
        Abstract method that must be implemented by subclasses. Called by ``parse_molecule()``
        to gather information for subclass-specific metrics
    parse_molecule(tags, record)
        Extract information from a set of sequencing reads that correspond to a molecule and store
        the data in the MetricAggregator class.
    finalize()
        Some metrics cannot be calculated until all the information for an entity has been
        aggregated, for example, the number of `fragments_per_molecule`. Finalize calculates all
        such higher-order metrics

    """

    def __init__(self):

        # type definitions
        Chromosome: int
        Strand: bool  # reverse = True, see pysam.AlignedSegment.is_reverse
        Position: int
        Fragment: Tuple[Chromosome, Position, Strand]

        # count information
        self.n_reads: int = 0
        self.noise_reads: int = 0  # long polymers, N-sequences; NotImplemented
        self._fragment_histogram: Counter[Fragment] = Counter()
        self._molecule_histogram: Counter[str] = Counter()

        # molecule information
        self._molecule_barcode_fraction_bases_above_30 = OnlineGaussianSufficientStatistic()
        self.perfect_molecule_barcodes = 0

        self._genomic_reads_fraction_bases_quality_above_30 = OnlineGaussianSufficientStatistic()
        self._genomic_read_quality = OnlineGaussianSufficientStatistic()

        # alignment location information
        self.reads_mapped_exonic = 0
        self.reads_mapped_intronic = 0
        self.reads_mapped_utr = 0

        # todo implement this once we have a gene model
        # self.reads_mapped_outside_window = 0  # reads should be within 1000 bases of UTR
        # self._read_distance_from_termination_site = OnlineGaussianSufficientStatistic()

        # alignment uniqueness information
        self.reads_mapped_uniquely = 0
        self.reads_mapped_multiple = 0
        self.duplicate_reads = 0

        # alignment splicing information
        self.spliced_reads = 0
        self.antisense_reads = 0
        self._plus_strand_reads = 0  # strand balance  # todo implement property here

        # higher-order methods, filled in by finalize() when all data is extracted
        self.molecule_barcode_fraction_bases_above_30_mean: float = None
        self.molecule_barcode_fraction_bases_above_30_variance: float = None
        self.genomic_reads_fraction_bases_quality_above_30_mean: float = None
        self.genomic_reads_fraction_bases_quality_above_30_variance: float = None
        self.genomic_read_quality_mean: float = None
        self.genomic_read_quality_variance: float = None
        self.n_molecules: float = None
        self.n_fragments: float = None
        self.reads_per_molecule: float = None
        self.reads_per_fragment: float = None
        self.fragments_per_molecule: float = None
        self.fragments_with_single_read_evidence: int = None
        self.molecules_with_single_read_evidence: int = None

    @staticmethod
    def _quality_string_to_numeric(quality_sequence: Iterable[str]) -> List[int]:
        """Convert an HTSlib ASCII quality string to an integer representation.

        Parameters
        ----------
        quality_sequence : Iterable[str]
            An iterable of Illumina base call qualities in ASCII encoding

        Returns
        -------
        numeric_qualities : List[int]
            A list of Illumina base call qualities converted to integers

        """
        return [ord(c) - 33 for c in quality_sequence]  # todo look up if this is accurate

    @staticmethod
    def _quality_above_threshold(threshold: int, quality_sequence: Sequence[int]) -> float:
        """Calculate the fraction of bases called with a quality above ``threshold``.

        Parameters
        ----------
        threshold: int
            The quality threshold
        quality_sequence: Sequence[int]
            A sequence of Illumina base qualities

        Returns
        -------
        fraction : float
            The fraction of bases in ``quality_sequence`` with quality greater than ``threshold``

        """
        return sum(1 for base in quality_sequence if base > threshold) / len(quality_sequence)

    def _is_noise(self, record: pysam.AlignedSegment) -> bool:
        return NotImplemented  # todo required because 10x measures this

    def parse_molecule(
            self, tags: Sequence[str], records: Iterable[pysam.AlignedSegment]) -> None:
        """Parse information from all records of a molecule.

        The parsed information is stored in the MetricAggregator in-place.

        Parameters
        ----------
        tags : Sequence[str]
            all the tags that define this molecule. one of {[CB, GE, UB], [GE, CB, UB]}
        records : Iterable[pysam.AlignedSegment]
            the sam records associated with the molecule

        """
        for record in records:

            # todo think about how I could use the duplicate tag to reduce computation; duplicates
            # should normally come in order in a sorted file

            # extract sub-class-specific information
            self.parse_extra_fields(tags=tags, record=record)

            self.n_reads += 1
            # self.noise_reads += self.is_noise(record)  # todo implement me

            # the tags passed to this function define a molecule, this increments the counter,
            # identifying a new molecule only if a new tag combination is observed
            self._molecule_histogram[tags] += 1

            self._molecule_barcode_fraction_bases_above_30.update(
                self._quality_above_threshold(
                    30, self._quality_string_to_numeric(record.get_tag(consts.QUALITY_MOLECULE_BARCODE_TAG_KEY))))

            # we should be tolerant and handle it if the pysam.AlignedSegment.get_tag 
            # cannot retrieve the data by a tag since it's not a fatal error
            try:
                self.perfect_molecule_barcodes += (
                    record.get_tag(consts.RAW_MOLECULE_BARCODE_TAG_KEY) == record.get_tag(consts.MOLECULE_BARCODE_TAG_KEY))
            except KeyError:
                # An error occurred while retrieving the data from the optional alighment section, which 
                # indicates that the read did not have a corrected UMI sequenct. In the future we would like to 
                # keep track of these reads.
                pass

            self._genomic_reads_fraction_bases_quality_above_30.update(
                self._quality_above_threshold(30, record.query_alignment_qualities))

            mean_alignment_quality: float = np.mean(record.query_alignment_qualities)
            self._genomic_read_quality.update(mean_alignment_quality)

            # the remaining portions deal with aligned reads, so if the read is not mapped, we are
            # done with it
            if record.is_unmapped:
                continue

            # get components that define a unique sequence fragment and increment the histogram
            position: int = record.pos
            strand: bool = record.is_reverse
            reference: int = record.reference_id
            self._fragment_histogram[reference, position, strand, tags] += 1

            alignment_location = record.get_tag(consts.ALIGNMENT_LOCATION_TAG_KEY)
            if alignment_location == consts.CODING_ALIGNMENT_LOCATION_TAG_VALUE:
                self.reads_mapped_exonic += 1
            elif alignment_location == consts.INTRONIC_ALIGNMENT_LOCATION_TAG_VALUE:
                self.reads_mapped_intronic += 1
            elif alignment_location == consts.UTR_ALIGNMENT_LOCATION_TAG_VALUE:
                self.reads_mapped_utr += 1

            # todo check if read maps outside window (needs gene model)
            # todo create distances from terminate side (needs gene model)

            # uniqueness
            number_mappings = record.get_tag(consts.NUMBER_OF_HITS_TAG_KEY)
            if number_mappings == 1:
                self.reads_mapped_uniquely += 1
            else:
                self.reads_mapped_multiple += 1  # todo without multi-mapping, this number is zero!

            if record.is_duplicate:
                self.duplicate_reads += 1

            # cigar N field (3) indicates a read is spliced if the value is non-zero
            cigar_stats, num_blocks = record.get_cigar_stats()
            if cigar_stats[3]:
                self.spliced_reads += 1

            # todo figure out antisense and make this notation clearer; info likely in dropseqtools
            self._plus_strand_reads += not record.is_reverse

    def parse_extra_fields(self, tags: Sequence[str], record: pysam.AlignedSegment) -> None:
        """Defined by subclasses to extract class-specific information from molecules"""
        raise NotImplementedError

    def finalize(self) -> None:
        """Calculate metrics that require information from all molecules of an entity

        ``finalize()`` replaces attributes in-place that were initialized by the constructor as
        ``None`` with a value calculated across all molecule data that has been aggregated.

        """

        self.molecule_barcode_fraction_bases_above_30_mean: float = \
            self._molecule_barcode_fraction_bases_above_30.mean

        self.molecule_barcode_fraction_bases_above_30_variance: float = \
            self._molecule_barcode_fraction_bases_above_30.calculate_variance()

        self.genomic_reads_fraction_bases_quality_above_30_mean: float = \
            self._genomic_reads_fraction_bases_quality_above_30.mean

        self.genomic_reads_fraction_bases_quality_above_30_variance: float = \
            self._genomic_reads_fraction_bases_quality_above_30.calculate_variance()

        self.genomic_read_quality_mean: float = self._genomic_read_quality.mean

        self.genomic_read_quality_variance: float = self._genomic_read_quality.calculate_variance()

        self.n_molecules: int = len(self._molecule_histogram.keys())

        self.n_fragments: int = len(self._fragment_histogram.keys())

        try:
            self.reads_per_molecule: float = self.n_reads / self.n_molecules
        except ZeroDivisionError:
            self.reads_per_molecule: float = float('nan')

        try:
            self.reads_per_fragment: float = self.n_reads / self.n_fragments
        except ZeroDivisionError:
            self.reads_per_fragment: float = float('nan')

        try:
            self.fragments_per_molecule: float = self.n_fragments / self.n_molecules
        except ZeroDivisionError:
            self.fragments_per_molecule: float = float('nan')

        self.fragments_with_single_read_evidence: int = \
            sum(1 for v in self._fragment_histogram.values() if v == 1)

        self.molecules_with_single_read_evidence: int = \
            sum(1 for v in self._molecule_histogram.values() if v == 1)


class CellMetrics(MetricAggregator):
    """Cell Metric Aggregator

    Aggregator that captures metric information about a cell by parsing all of the molecules in
    an experiment that were annotated with a specific cell barcode, as recorded in the ``CB`` tag.

    Attributes
    ----------
    perfect_cell_barcodes : int
        The number of reads whose cell barcodes contain no errors (tag ``CB`` == ``CR``)
    reads_mapped_intergenic : int
        The number of reads mapped to an intergenic region for this cell
    reads_mapped_too_many_loci : int
        The number of reads that were mapped to too many loci across the genome and as a
        consequence, are reported unmapped by the aligner
    cell_barcode_fraction_bases_above_30_variance : float
        The variance of the fraction of Illumina base calls for the cell barcode sequence that
        are greater than 30, across molecules
    cell_barcode_fraction_bases_above_30_mean : float
        The average fraction of Illumina base calls for the cell barcode sequence that
        are greater than 30, across molecules
    n_genes : int
        The number of genes detected by this cell
    genes_detected_multiple_observations : int
        The number of genes that are observed by more than one read in this cell

    """

    extra_docs = """
    Examples
    --------
    # todo implement me

    See Also
    --------
    GeneMetrics

    """

    __doc__ += MetricAggregator.__doc__ + extra_docs

    def __init__(self):
        super().__init__()

        # barcode quality data
        self._cell_barcode_fraction_bases_above_30 = OnlineGaussianSufficientStatistic()
        self.perfect_cell_barcodes = 0  # inv: fraction cells with errors

        # track non-transcriptomic reads
        self.reads_mapped_intergenic = 0
        self.reads_unmapped = 0
        self.reads_mapped_too_many_loci = 0

        self._genes_histogram = Counter()

        # todo think about whether we can build molecule models that map to things that aren't genes
        # i.e. to integentic regions or intronic regions. This could be a part of multi-mapping
        # self.molecules_mapped_intergenic = 0

        self.cell_barcode_fraction_bases_above_30_variance: float = None
        self.cell_barcode_fraction_bases_above_30_mean: float = None
        self.n_genes: int = None
        self.genes_detected_multiple_observations: int = None

    def finalize(self):
        super().finalize()

        self.cell_barcode_fraction_bases_above_30_mean: float = \
            self._cell_barcode_fraction_bases_above_30.mean

        self.cell_barcode_fraction_bases_above_30_variance: float = \
            self._cell_barcode_fraction_bases_above_30.calculate_variance()

        self.n_genes: int = \
            len(self._genes_histogram.keys())

        self.genes_detected_multiple_observations: int = \
            sum(1 for v in self._genes_histogram.values() if v > 1)

    def parse_extra_fields(self, tags: Sequence[str], record: pysam.AlignedSegment) -> None:
        """Parses a record to extract gene-specific information

        Gene-specific metric data is stored in-place in the MetricAggregator

        Parameters
        ----------
        tags : Sequence[str]
            The GE, UB and CB tags that define this molecule
        record : pysam.AlignedSegment
            SAM record to be parsed

        """
        self._cell_barcode_fraction_bases_above_30.update(
            self._quality_above_threshold(
                30, self._quality_string_to_numeric(record.get_tag(consts.QUALITY_CELL_BARCODE_TAG_KEY))))

        # Exclude reads that do not have a CB tag from the perfect_cell_barcodes count
        if record.has_tag(consts.CELL_BARCODE_TAG_KEY):
            raw_cell_barcode_tag = record.get_tag(consts.RAW_CELL_BARCODE_TAG_KEY)
            cell_barcode_tag = record.get_tag(consts.CELL_BARCODE_TAG_KEY)
            self.perfect_cell_barcodes += (raw_cell_barcode_tag == cell_barcode_tag)

        try:
            alignment_location = record.get_tag(consts.ALIGNMENT_LOCATION_TAG_KEY)
            if alignment_location == consts.INTERGENIC_ALIGNMENT_LOCATION_TAG_VALUE:
                self.reads_mapped_intergenic += 1
        except KeyError:
            self.reads_unmapped += 1

        # todo track reads_mapped_too_many_loci after multi-alignment is done

        self._genes_histogram[tags[2]] += 1  # note that no gene == None


class GeneMetrics(MetricAggregator):
    """Gene Metric Aggregator

    Aggregator that captures metric information about a gene by parsing all of the molecules in
    an experiment that were annotated with a specific gene ID, as recorded in the ``GE`` tag.

    Attributes
    ----------
    number_cells_detected_multiple : int
        The number of cells which observe more than one read of this gene
    number_cells_expressing : int
        The number of cells that detect this gene

    """

    extra_docs = """
    Examples
    --------
    # todo implement me

    See Also
    --------
    CellMetrics

    """

    __doc__ += MetricAggregator.__doc__ + extra_docs

    def __init__(self):
        super().__init__()

        self._cells_histogram = Counter()
        # todo we don't tag exon right now. Not sure if we want to or not
        # self._exon_histogram = Counter()

        self.number_cells_detected_multiple: int = None
        self.number_cells_expressing: int = None

    def finalize(self):
        super().finalize()

        self.number_cells_expressing: int = len(self._cells_histogram.keys())

        self.number_cells_detected_multiple: int = \
            sum(1 for c in self._cells_histogram.values() if c > 1)

    def parse_extra_fields(self, tags: Sequence[str], record: pysam.AlignedSegment) -> None:
        """Parses a record to extract cell-specific information

        Cell-specific metric data is stored in-place in the MetricAggregator

        Parameters
        ----------
        tags : Sequence[str]
            The CB, UB and GE tags that define this molecule
        record : pysam.AlignedSegment
            SAM record to be parsed

        """
        self._cells_histogram[tags[1]] += 1
