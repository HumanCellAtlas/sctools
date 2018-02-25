from typing import TextIO
from contextlib import closing
import pysam
from sctools.bam import iter_cell_barcodes, iter_genes, iter_molecule_barcodes
from sctools.metrics.aggregator import CellBarcodeMetrics, GeneMetrics


class MetricGatherer:

    def __init__(self, bam_file: str, output_stem: str):
        self._bam_file = bam_file
        self._output_stem = output_stem

    @property
    def bam_file(self) -> str:
        return self._bam_file

    def generate_output_file(self) -> TextIO:
        if self._output_stem.endswith('.csv'):
            return open(self._output_stem, 'w')
        else:
            return open(self._output_stem + '.csv', 'w')

    def extract_metrics(self):
        raise NotImplementedError


class GatherCellMetrics(MetricGatherer):

    # requires that bam file is sorted by GE, UB, CB tags

    def extract_metrics(self) -> None:

        # open the files
        with pysam.AlignmentFile(self.bam_file, 'rb') as bam_iterator, \
                closing(self.generate_output_file()) as cell_metrics_output:

            # write the header
            CellBarcodeMetrics().write_csv_header(cell_metrics_output)

            # break up the bam file into sub-iterators over cells
            for cell_iterator, cell_tag in iter_cell_barcodes(bam_iterator=bam_iterator):
                metric_aggregator = CellBarcodeMetrics()

                # break up cell barcodes by molecule barcodes
                for molecule_iterator, molecule_tag in iter_molecule_barcodes(
                        bam_iterator=cell_iterator):

                    # todo there are a bunch of reads that don't have genes that this misses
                    # break up molecule barcodes by gene ids
                    for gene_iterator, gene_tag in iter_genes(bam_iterator=molecule_iterator):

                        # process the data
                        metric_aggregator.parse_molecule(
                            tags=(cell_tag, molecule_tag, gene_tag),
                            records=gene_iterator
                        )

                # write a record for each cell
                metric_aggregator.write_csv_record(open_fid=cell_metrics_output, index=cell_tag)


class GatherGeneMetrics(MetricGatherer):

    # assumes file is sorted by UB, CB, GE tag

    def extract_metrics(self):
        # open the files
        with pysam.AlignmentFile(self.bam_file, 'rb') as bam_iterator, \
                closing(self.generate_output_file()) as gene_metrics_output:

            # write the header
            GeneMetrics().write_csv_header(gene_metrics_output)

            # break up the bam file into sub-iterators over cells
            for gene_iterator, gene_tag in iter_genes(bam_iterator=bam_iterator):
                metric_aggregator = GeneMetrics()

                # break up cell barcodes by molecule barcodes
                for cell_iterator, cell_tag in iter_cell_barcodes(bam_iterator=gene_iterator):

                    # break up molecule barcodes by gene ids
                    for molecule_iterator, molecule_tag in iter_molecule_barcodes(
                            bam_iterator=cell_iterator):

                        # process the data
                        metric_aggregator.parse_molecule(
                            tags=(gene_tag, cell_tag, molecule_tag),
                            records=molecule_iterator
                        )

                # write a record for each cell
                metric_aggregator.finalize()
                metric_aggregator.write_csv_record(open_fid=gene_metrics_output, index=gene_tag)
