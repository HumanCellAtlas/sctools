from typing import Iterator, Iterable, TextIO
from collections import Counter, OrderedDict
import pysam


class MetricAggregator:

    _data: OrderedDict = NotImplemented

    @property
    def data(self) -> OrderedDict:
        return self._data

    def write_csv_header(self, open_fid):
        open_fid.write(',' + ','.join(self.data.keys()) + '\n')

    def write_csv_record(self, open_fid):
        open_fid.write(','.join(self.data.values()) + '\n')

    def parse_records(self, records: Iterable[pysam.AlignedSegment]) -> None:
        """
        this function extracts information from all the records of a molecule and adds them to the
        parent aggregator
        """
        raise NotImplementedError


class CellBarcodeMetrics(MetricAggregator):

    def __init__(self):
        self._data = OrderedDict(
            # count information
            number_reads=None,
            number_molecules=None,
            # dependent metric: reads / molecule

            # barcode quality data
            number_noise_reads=None,  # long polymers, N-sequences
            number_cell_barcode_bases_above_30=None,
            number_perfect_cell_barcodes=None,  # inv: fraction cells with errors
            average_number_molecule_barcode_bases_above_30=None,
            number_perfect_molecule_barcodes=None,  # inv: fraction molecules with errors
            average_number_genomic_reads_quality_above_30=None,

            # alignment locations
            number_reads_mapped_to_UTR=None,
            number_reads_mapped_to_transcriptome=None,
            number_reads_mapped_intergenic=None,
            number_reads_mapped_intronic=None,
            number_reads_mapped_uniquely=None,
            number_reads_mapped_too_many_loci=None,
            number_reads_mapped_multiple=None,
            number_reads_mapped_outside_window=None,  # reads should be within 1000 bases of UTR
            number_duplicate_reads=None,
            number_genes_detected=None,

            plus_strand_reads=None,  # strand balance

            number_spliced_reads=None,
            number_antisense_reads=None,

            average_number_fragments_per_molecule=None,
            number_molecules_with_single_read_evidence=None,

            number_molecules_mapped_to_UTR=None,
            number_molecules_mapped_to_transcriptome=None,
            number_molecules_mapped_intergenic=None,
            number_molecules_mapped_intronic=None,

            # higher-order statistics:
            yield_=None,
        )

    def parse_records(self, records: Iterable[pysam.AlignedSegment]) -> None:
        """
        this function extracts information from all the records of a molecule and adds them to the
        parent aggregator
        """
        raise NotImplementedError  # implement this


class Gene(MetricAggregator):

    def __init__(self):
        self._data = OrderedDict(
            # read metrics
            number_reads_mapped_to_UTR=None,
            number_reads_mapped_to_transcriptome=None,
            number_reads_mapped_intronic=None,
            number_reads_mapped_uniquely=None,
            number_reads_mapped_multiple=None,
            number_reads_mapped_outside_window=None,  # reads should be within 1000 bases of UTR
            number_duplicate_reads=None,
            average_distance_from_termination_site=None,

            # strand and splicing
            plus_strand_reads=None,  # strand balance
            number_spliced_reads=None,
            number_antisense_reads=None,

            # molecule confidence measurements
            average_number_fragments_per_molecule=None,
            number_molecules_with_single_read_evidence=None,

            # molecule metrics
            number_molecules_mapped_to_UTR=None,
            number_molecules_mapped_to_transcriptome=None,
            number_molecules_mapped_intronic=None,
            number_distinct_mapped_UMIs=None,

            # higher-order statistics:
            yield_=None,
        )

        exon_histogram = Counter()

    def parse_records(self, record: Iterable[pysam.AlignedSegment]) -> None:
        """
        This function extracts information from all the records of a gene and adds them to the
        parent aggregator
        """
        raise NotImplementedError  # implement this


class MetricGatherer:

    def __init__(self, bam_file: str):
        self._bam_file = bam_file

    @staticmethod
    def group_by_tag(self, iterable: Iterable) -> Iterable:
        pass  # function to group an iterable by a tag and produce another iterable

    def generate_output_file(self) -> TextIO:
        pass  # some kind of CSV file, for portability

    def parse_bam(self):
        raise NotImplementedError


class GatherCellBarcodeMetrics(MetricGatherer):

    # assumes file is sorted by CB then UB tags

    def iter_cells(self, record_iterator: Iterator) -> Iterator:
        """
        function to iterate over all the cells of a bam file
        """
        raise NotImplementedError

    def iter_molecules(self, cell_record_iterator: Iterator) -> Iterator:
        """
        function to iterate over all the molecules of a cell
        """
        raise NotImplementedError

    def parse_bam(self):

        # wasn't sure where to put histograms; not sure if I can abstract much more here

        # histogram supports: cell entropy, number cell barcodes, valid barcodes
        cell_barcode_histogram = Counter()

        # histogram supports: umi entropy, priors for downstream analysis
        molecule_barcode_histogram = Counter()

        # initialize the aggregator
        cell_barcode_open_fid = self.generate_output_file()

        for cell in self.iter_cells(self._bam_file):
            aggregator = CellBarcodeMetrics()
            for molecule in self.iter_molecules(cell):
                aggregator.parse_records(molecule)
            aggregator.write_csv_record(cell_barcode_open_fid)

        # todo post-process the histograms for the metrics they contain


class GatherGeneMetrics(MetricGatherer):

    # assumes file is sorted by GE tag

    def iter_genes(self, bam_file) -> Iterator:
        raise NotImplementedError  # get all records for a gene

    def parse_bam(self):

        gene_open_fid = self.generate_output_file()

        for gene_records in self.iter_genes(self._bam_file):
            gene = Gene()
            gene.parse_records(gene_records)
            gene.write_csv_record(gene_open_fid)


class MergeMetrics:

    def __init__(self, metric_files: Iterable[str]):
        self._metric_files = metric_files

    def merge_function(self):
        raise NotImplementedError  # merge the metrics


class MergeCellMetrics(MergeMetrics):
    pass  # cell metrics can be concatenated


class MergeGeneMetrics(MergeMetrics):
    pass  # gene metrics must be averaged and merged

