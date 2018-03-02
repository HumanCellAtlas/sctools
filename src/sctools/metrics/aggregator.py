from typing import Iterable, Tuple, Counter, List, Sequence
from collections import Counter
import pysam
from sctools.stats import OnlineGaussianSufficientStatistic
import numpy as np

# note that this entire module can be rewritten with dataclass when python 3.7 stabilizes
# see https://www.python.org/dev/peps/pep-0557/
# the utility of this pep has driven some of the implementation details here


class SequenceMetricAggregator:

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

    def finalize(self):

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

    @staticmethod
    def quality_string_to_numeric(quality_sequence: Iterable[str]) -> List[int]:
        return [ord(c) - 33 for c in quality_sequence]  # todo look up if this is accurate

    @staticmethod
    def quality_above_threshold(threshold: int, quality_sequence: Sequence[int]) -> float:
        return sum(1 for base in quality_sequence if base > threshold) / len(quality_sequence)

    # def is_noise(self, record: pysam.AlignedSegment) -> bool:
    #     """returns True if the record is noise (polymeric or lots of Ns, read 10x code)"""
    #     raise NotImplementedError

    @staticmethod
    def is_duplicate(record: pysam.AlignedSegment) -> bool:
        return record.is_duplicate

    def parse_molecule(
            self, tags: Sequence[str], records: Iterable[pysam.AlignedSegment]) -> None:
        """Parse records of a molecule, storing the resulting data

        :param tags: all the tags that define this molecule
        :param records: the records associated with the molecule
        """
        for record in records:

            # todo think about how I could use the duplicate tag to simplify this class

            # call the parse_extra_fields function to extract subclass-specific information
            # todo might need a callback if this function has knowledge about termination conditions
            self.parse_extra_fields(tags=tags, record=record)

            self.n_reads += 1
            # self.noise_reads += self.is_noise(record)

            # get components that define a unique sequence fragment
            position: int = record.pos
            strand: bool = record.is_reverse
            reference: int = record.reference_id
            self._fragment_histogram[reference, position, strand, tags] += 1

            self._molecule_histogram[tags] += 1  # might also need gene here (or cell?)

            self._molecule_barcode_fraction_bases_above_30.update(
                self.quality_above_threshold(
                    30, self.quality_string_to_numeric(record.get_tag('UY'))
                )
            )

            self.perfect_molecule_barcodes += record.get_tag('UR') == record.get_tag('UB')

            self._genomic_reads_fraction_bases_quality_above_30.update(
                self.quality_above_threshold(30, record.query_alignment_qualities)
            )
            self._genomic_read_quality.update(
                np.mean(record.query_qualities)
            )

            # alignment location information
            alignment_location = record.get_tag('XF')
            if alignment_location == 'CODING':
                self.reads_mapped_exonic += 1
            elif alignment_location == 'INTRONIC':
                self.reads_mapped_intronic += 1
            elif alignment_location == 'UTR':
                self.reads_mapped_utr += 1

            # todo check if read maps outside window (need gene model)
            # todo create distances from terminate side (need gene model)

            # uniqueness
            number_mappings = record.get_tag('NH')
            if number_mappings == 1:
                self.reads_mapped_uniquely += 1
            else:
                self.reads_mapped_multiple += 1  # todo without MM, this number is zero!

            if self.is_duplicate(record):
                self.duplicate_reads += 1

            # splicing
            cigar_stats, num_blocks = record.get_cigar_stats()
            if cigar_stats[3]:
                self.spliced_reads += 1
            # todo figure out antisense -- this is probably dropseqtools somewhere
            self._plus_strand_reads += not record.is_reverse

    def parse_extra_fields(self, tags: Sequence[str], record: pysam.AlignedSegment) -> None:
        """
        this function extracts additional metrics specific to the parent extractor and adds them
        to the parent.
        """
        raise NotImplementedError


class CellBarcodeMetrics(SequenceMetricAggregator):

    def __init__(self):
        super().__init__()

        # barcode quality data
        self._cell_barcode_fraction_bases_above_30 = OnlineGaussianSufficientStatistic()
        self.perfect_cell_barcodes = 0  # inv: fraction cells with errors

        # track non-transcriptomic reads
        self.reads_mapped_intergenic = 0
        self.reads_mapped_too_many_loci = 0

        self._genes_histogram = Counter()

        # todo think about whether we can build molecule models that map to things that aren't genes
        # this could be a part of multi-mapping
        # self.molecules_mapped_intergenic = 0

        self.cell_barcode_fraction_bases_above_30_variance: float = None
        self.cell_barcode_fraction_bases_above_30_mean: float = None
        self.genes_detected: int = None
        self.genes_detected_multiple_observations: int = None

    def finalize(self):
        super().finalize()

        self.cell_barcode_fraction_bases_above_30_mean: float = \
            self._cell_barcode_fraction_bases_above_30.mean

        self.cell_barcode_fraction_bases_above_30_variance: float = \
            self._cell_barcode_fraction_bases_above_30.calculate_variance()

        self.genes_detected: int = \
            len(self._genes_histogram.keys())

        self.genes_detected_multiple_observations: int = \
            sum(1 for v in self._genes_histogram.values() if v > 1)

    def parse_extra_fields(self, tags: Sequence[str], record: pysam.AlignedSegment) -> None:
        """
        this function extracts metric information specific to this subclass and adds them to the
        parent aggregator. It is called by parse_molecule on each record.

        if this function returns false, the full suite is not run for this record
        """
        self._cell_barcode_fraction_bases_above_30.update(
            self.quality_above_threshold(
                30, self.quality_string_to_numeric(record.get_tag('CY'))
            )
        )

        self.perfect_cell_barcodes += record.get_tag('UR') == record.get_tag('UB')
        alignment_location = record.get_tag('XF')
        if alignment_location == 'INTERGENIC':
            self.reads_mapped_intergenic += 1

        # todo track reads_mapped_too_many_loci after multi-alignment is done

        self._genes_histogram[tags[0]] += 1


class GeneMetrics(SequenceMetricAggregator):

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
        """
        this function extracts metric information specific to this subclass and adds them to the
        parent aggregator. It is called by parse_molecule on each record
        """
        self._cells_histogram[tags[1]] += 1
