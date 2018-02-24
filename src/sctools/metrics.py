from typing import Iterator, Iterable, TextIO, Generator, List, Tuple
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

    def parse_records(self, records: Iterable[pysam.AlignedSegment], tag: str) -> None:
        """
        this function extracts information from all the records of a given tag or tag combination
        and adds them to the parent aggregator
        """
        raise NotImplementedError


class CellBarcodeMetrics(MetricAggregator):

    def __init__(self):
        self._data = OrderedDict(
            # count information
            reads=0,
            molecules=0,
            # dependent metric: reads / molecule

            # barcode quality data
            noise_reads=0,  # long polymers, N-sequences
            cell_barcode_bases_above_30=0,
            perfect_cell_barcodes=0,  # inv: fraction cells with errors
            average_molecule_barcode_bases_above_30=0,
            perfect_molecule_barcodes=0,  # inv: fraction molecules with errors
            average_genomic_reads_quality_above_30=0,

            # alignment locations
            reads_mapped_to_UTR=0,
            reads_mapped_to_transcriptome=0,
            reads_mapped_intergenic=0,
            reads_mapped_intronic=0,
            reads_mapped_uniquely=0,
            reads_mapped_too_many_loci=0,
            reads_mapped_multiple=0,
            reads_mapped_outside_window=0,  # reads should be within 1000 bases of UTR
            duplicate_reads=0,
            genes_detected=0,

            plus_strand_reads=0,  # strand balance

            spliced_reads=0,
            antisense_reads=0,

            average_fragments_per_molecule=0,
            molecules_with_single_read_evidence=0,

            molecules_mapped_to_UTR=0,
            molecules_mapped_to_transcriptome=0,
            molecules_mapped_intergenic=0,
            molecules_mapped_intronic=0,
        )

    def calculate_higher_order_metrics(self):
        """calculate any metrics that depend on aggregated data"""
        pass  # yield

    def parse_records(self, records: Iterable[pysam.AlignedSegment], tag: str) -> None:
        """
        this function extracts information from all the records of a molecule and adds them to the
        parent aggregator
        """
        raise NotImplementedError  # implement this


class Gene(MetricAggregator):

    def __init__(self):
        self._data = OrderedDict(
            # read metrics
            reads_unmapped=0,
            reads_mapped_to_UTR=0,
            reads_mapped_to_transcriptome=0,
            reads_mapped_intronic=0,
            reads_mapped_uniquely=0,
            reads_mapped_multiple=0,
            reads_mapped_outside_window=0,  # reads should be within 1000 bases of UTR
            duplicate_reads=0,
            average_distance_from_termination_site=0,

            # strand and splicing
            plus_strand_reads=0,  # strand balance
            spliced_reads=0,
            antisense_reads=0,

            # molecule confidence measurements
            average_fragments_per_molecule=0,
            molecules_with_single_read_evidence=0,
            average_read_quality=0,

            # molecule metrics
            molecules_mapped_to_UTR=0,
            molecules_mapped_to_transcriptome=0,
            molecules_mapped_intronic=0,
            distinct_mapped_UMIs=0,
        )

        exon_histogram = Counter()

    def calculate_higher_order_metrics(self):
        """calculate any metrics that depend on aggregated data"""
        pass  # yield

    def parse_records(self, records: Iterable[pysam.AlignedSegment], tag: str) -> None:
        """
        This function extracts information from all the records of a gene and adds them to the
        parent aggregator
        """
        for record in records:
            # check if/where the record maps
            if record.is_unmapped:
                self.data['reads_unmapped'] += 1
                return
            else:
                record.get_tag('XF')  # stopped here


class MetricGatherer:

    def __init__(self, bam_file: str, output_stem: str):
        self._bam_file = bam_file
        self._output_stem = output_stem

    @property
    def bam_file(self) -> str:
        return self._bam_file

    @staticmethod
    def iter_tag_groups(tag: str, bam_iterator: Iterator[pysam.AlignedSegment]) -> Generator[List]:

        # get first read and tag set
        reads = [next(bam_iterator)]
        current_tag = reads[0].get_tag(tag)

        for alignment in bam_iterator:
            next_tag = alignment.get_tag(tag)
            if next_tag == current_tag:
                reads.append(alignment)
            else:
                yield iter(reads), current_tag
                reads = [alignment]
                current_tag = next_tag
        yield iter(reads), current_tag

    def generate_output_file(self) -> TextIO:
        return open(self._output_stem, 'w')

    def parse_bam(self):
        raise NotImplementedError


class GatherCellMetrics(MetricGatherer):

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

    def parse_bam(self) -> None:

        # wasn't sure where to put histograms; not sure if I can abstract much more here

        # histogram supports: cell entropy, number cell barcodes, valid barcodes
        cell_barcode_histogram = Counter()

        # histogram supports: umi entropy, priors for downstream analysis
        molecule_barcode_histogram = Counter()

        # initialize the aggregator
        cell_metrics = self.generate_output_file()

        for cell, cell_tag in self.iter_tag_groups('CB', self._bam_file):
            aggregator = CellBarcodeMetrics()
            cell_barcode_histogram[cell_tag] += 1
            for molecule, molecule_tag in self.iter_tag_groups('UB', cell):
                aggregator.parse_records(molecule)
                molecule_barcode_histogram[molecule_tag] += 1
            aggregator.write_csv_record(cell_metrics)

        # todo post-process the histograms for the metrics they contain


class GatherGeneMetrics(MetricGatherer):

    # assumes file is sorted by GE tag

    def parse_bam(self):

        gene_open_fid = self.generate_output_file()

        for gene_records in self.iter_tag_groups('GE', self._bam_file):
            gene = Gene()
            gene.parse_records(gene_records)
            gene.write_csv_record(gene_open_fid)


class MergeMetrics:

    def __init__(self, metric_files: Iterable[str], output_file: str):
        self._metric_files = metric_files
        self._output_file = output_file

    def merge_function(self):
        raise NotImplementedError  # merge the metrics


class MergeCellMetrics(MergeMetrics):
    pass  # cell metrics can be concatenated


class MergeGeneMetrics(MergeMetrics):
    pass  # gene metrics must be averaged and merged

