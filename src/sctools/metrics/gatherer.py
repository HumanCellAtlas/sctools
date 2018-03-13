from contextlib import closing
import pysam
from sctools.bam import iter_cell_barcodes, iter_genes, iter_molecule_barcodes
from sctools.metrics.aggregator import CellBarcodeMetrics, GeneMetrics
from sctools.metrics.writer import MetricCSVWriter


class MetricGatherer:

    def __init__(self, bam_file: str, output_stem: str):
        """
        MetricGatherer defines an extract_metrics() function which will parse a file molecule-by-
        molecule and write the resulting information to an output file.

        Because molecules tend to have relatively small numbers of reads, the memory footprint of
        this method is typically small (tens of megabytes)

        :param bam_file: the bam file over which to compile metrics
        :param output_stem: the file stem for the (csv formatted) metric outputs
        """
        self._bam_file = bam_file
        self._output_stem = output_stem

    @property
    def bam_file(self) -> str:
        return self._bam_file

    def extract_metrics(self) -> None:
        """extract metrics from the provided bam file and write the results to csv."""
        raise NotImplementedError


class GatherCellMetrics(MetricGatherer):

    # assumes bam file is sorted by GE, UB, CB tags

    def extract_metrics(self, mode: str='rb') -> None:

        # open the files
        with pysam.AlignmentFile(self.bam_file, mode=mode) as bam_iterator, \
                closing(MetricCSVWriter(self._output_stem)) as cell_metrics_output:

            # write the header
            cell_metrics_output.write_header(vars(CellBarcodeMetrics()))

            # break up the bam file into sub-iterators over cell barcodes
            for cell_iterator, cell_tag in iter_cell_barcodes(bam_iterator=bam_iterator):
                metric_aggregator = CellBarcodeMetrics()

                # break up cell barcodes by molecule barcodes
                for molecule_iterator, molecule_tag in iter_molecule_barcodes(
                        bam_iterator=cell_iterator):

                    # break up molecule barcodes by gene ids
                    for gene_iterator, gene_tag in iter_genes(bam_iterator=molecule_iterator):

                        # process the data
                        metric_aggregator.parse_molecule(
                            tags=(cell_tag, molecule_tag, gene_tag),
                            records=gene_iterator
                        )

                # write a record for each cell
                metric_aggregator.finalize()
                cell_metrics_output.write(cell_tag, vars(metric_aggregator))


class GatherGeneMetrics(MetricGatherer):

    # assumes bam file is sorted by UB, CB, GE tag

    def extract_metrics(self, mode: str='rb'):
        # open the files
        with pysam.AlignmentFile(self.bam_file, mode=mode) as bam_iterator, \
                closing(MetricCSVWriter(self._output_stem)) as gene_metrics_output:

            # write the header
            gene_metrics_output.write_header(vars(GeneMetrics()))

            # break up the bam file into sub-iterators over gene ids
            for gene_iterator, gene_tag in iter_genes(bam_iterator=bam_iterator):
                metric_aggregator = GeneMetrics()

                # break up gene ids by cell barcodes
                for cell_iterator, cell_tag in iter_cell_barcodes(bam_iterator=gene_iterator):

                    # break up cell barcodes by molecular barcodes
                    for molecule_iterator, molecule_tag in iter_molecule_barcodes(
                            bam_iterator=cell_iterator):

                        # process the data
                        metric_aggregator.parse_molecule(
                            tags=(gene_tag, cell_tag, molecule_tag),
                            records=molecule_iterator
                        )

                # write a record for each gene id
                metric_aggregator.finalize()
                gene_metrics_output.write(gene_tag, vars(metric_aggregator))
