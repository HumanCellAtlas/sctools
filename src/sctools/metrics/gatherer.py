"""
Sequence Metric Gatherers
=========================

..currentmodule:: sctools.metrics

This module defines classes to gather metrics across the cells or genes of an experiment and write
them to gzip-compressed csv files

Classes
-------

.. autosummary::
   :toctree: generated/

   MetricGatherer               Gatherer Base Class
   GatherCellMetrics            Class to gather metrics on all cells in an experiment
   GatherGeneMetrics            Class to gather metrics on all genes in an experiment

See Also
--------
sctools.metrics.aggregator
sctools.metrics.merge
sctools.metrics.writer

"""

from contextlib import closing

import pysam
from typing import Set

from sctools.bam import iter_cell_barcodes, iter_genes, iter_molecule_barcodes
from sctools.metrics.aggregator import CellMetrics, GeneMetrics
from sctools.metrics.writer import MetricCSVWriter
from sctools.tsv import iter_tag_groups_from_tsv


class MetricGatherer:
    """Gathers Metrics from an experiment

    Because molecules tend to have relatively small numbers of reads, the memory footprint of
    this method is typically small (tens of megabytes).

    Parameters
    ----------
    bam_file : str
        the bam file containing the reads that metrics should be calculated from. Can be a chunk
        of cells or an entire experiment
    output_stem : str
        the file stem for the gzipped csv output

    Methods
    -------
    extract_metrics
        extracts metrics from ``bam_file`` and writes them to output_stem.csv.gz

    """

    def __init__(
        self,
        bam_file: str,
        output_stem: str,
        mitochondrial_gene_ids: Set[str] = set(),
        compress: bool = True,
    ):
        self._bam_file = bam_file
        self._output_stem = output_stem
        self._compress = compress
        self._mitochondrial_gene_ids = mitochondrial_gene_ids

    @property
    def bam_file(self) -> str:
        """the bam file that metrics are generated from"""
        return self._bam_file

    def extract_metrics(self, mode="rb") -> None:
        """extract metrics from the provided bam file and write the results to csv.

        Parameters
        ----------
        mode : {'r', 'rb'}, default 'rb'
            the open mode for pysam.AlignmentFile. 'r' indicates the input is a sam file, and 'rb'
            indicates a bam file.

        """
        raise NotImplementedError

    def compute_metrics(self, entitywise_records, metric_aggregator):
        # create a metrics aggregator for the current cell
        # process the records of the new cell
        prev_second_third_tags = (None, None)
        prev_second_tag = None
        prev_third_tag = None
        cellwise_data_field_records = []

        for _entitywise_record in entitywise_records:
            entitywise_record = [x.strip() for x in _entitywise_record.split("\t")]
            first_tag, second_tag, third_tag = (
                entitywise_record[0],
                entitywise_record[1],
                entitywise_record[2],
            )

            second_third_tags = (second_tag, third_tag)
            # not the first record and either molecule tar or gene tag is different from previous
            if (
                prev_second_third_tags != second_third_tags
                and cellwise_data_field_records != []
            ):

                # compute metrics for the molecule and gene tags
                metric_aggregator.parse_molecule_fast(
                    tags=(first_tag, prev_second_tag, prev_third_tag),
                    records=cellwise_data_field_records,
                )
                #                             print(cell_tag, molecule_tag, gene_tag)
                # create for the next combination of molecule and gene_tag records
                cellwise_data_field_records = []

            # otherwise continue with the data
            cellwise_data_field_records.append(
                entitywise_record[3:]
            )  # the first two are molecute and gene_tags
            prev_second_third_tags = second_third_tags
            prev_second_tag = second_tag
            prev_third_tag = third_tag

        #  Now process the last batch
        if cellwise_data_field_records:
            metric_aggregator.parse_molecule_fast(
                tags=(first_tag, prev_second_tag, prev_third_tag),
                records=cellwise_data_field_records,
            )
        cellwise_data_field_records = []
        return first_tag

class GatherCellMetrics(MetricGatherer):

    extra_docs = """
    Notes
    -----
    ``bam_file`` must be sorted by gene (``GE``), molecule (``UB``), and cell (``CB``), where gene
    varies fastest.

    Examples
    --------
    >>> from sctools.metrics.gatherer import GatherCellMetrics
    >>> import os, tempfile

    >>> # example data
    >>> bam_file = os.path.abspath(__file__) + '../test/data/test.bam'
    >>> temp_dir = tempfile.mkdtemp()
    >>> g = GatherCellMetrics(bam_file=bam_file, output_stem=temp_dir + 'test', compress=True)
    >>> g.extract_metrics()

    See Also
    --------
    GatherGeneMetrics

    """

    __doc__ += extra_docs

    def extract_metrics(self, mode: str = "rb") -> None:
        """Extract cell metrics from self.bam_file

        Parameters
        ----------
        mode : str, optional
            Open mode for self.bam. 'r' -> sam, 'rb' -> bam (default = 'rb').

        """
        # open the files
        with pysam.AlignmentFile(self.bam_file, mode=mode) as bam_iterator, closing(
            MetricCSVWriter(self._output_stem, self._compress)
        ) as cell_metrics_output:

            # write the header
            cell_metrics_output.write_header(vars(CellMetrics()))

            # break up the bam file into sub-iterators over cell barcodes
            for cell_iterator, cell_tag in iter_cell_barcodes(
                bam_iterator=bam_iterator
            ):
                metric_aggregator = CellMetrics()

                # break up cell barcodes by molecule barcodes
                for molecule_iterator, molecule_tag in iter_molecule_barcodes(
                    bam_iterator=cell_iterator
                ):

                    # break up molecule barcodes by gene ids
                    for gene_iterator, gene_tag in iter_genes(
                        bam_iterator=molecule_iterator
                    ):

                        # process the data
                        metric_aggregator.parse_molecule(
                            tags=(cell_tag, molecule_tag, gene_tag),
                            records=gene_iterator,
                        )

                # write a record for each cell
                metric_aggregator.finalize(
                    mitochondrial_genes=self._mitochondrial_gene_ids
                )
                cell_metrics_output.write(cell_tag, vars(metric_aggregator))

class GatherCellMetricsFast(MetricGatherer):

    extra_docs = """
    Notes
    -----
    ``bam_file`` must be sorted by gene (``GE``), molecule (``UB``), and cell (``CB``), where gene
    varies fastest.

    Examples
    --------
    >>> from sctools.metrics.gatherer import GatherCellMetrics
    >>> import os, tempfile

    >>> # example data
    >>> bam_file = os.path.abspath(__file__) + '../test/data/test.bam'
    >>> temp_dir = tempfile.mkdtemp()
    >>> g = GatherCellMetrics(bam_file=bam_file, output_stem=temp_dir + 'test', compress=True)
    >>> g.extract_metrics()

    See Also
    --------
    GatherGeneMetrics

    """

    __doc__ += extra_docs

    def extract_metrics(self, mode: str = "r") -> None:
        """Extract cell metrics from self.bam_file

        Parameters
        ----------
        mode : str, optional
            Open mode for self.bam. 'r' -> sam, 'rb' -> bam (default = 'rb').

        """

        # open the files
        with gzip.open(self.bam_file, mode=mode) if self.bam_file.endswisth('.gz') else \
            open(self.bam_file, mode=mode)  as tsv_reader, closing(
            MetricCSVWriter(self._output_stem, self._compress)
        ) as cell_metrics_output:
            # write the header
            cell_metrics_output.write_header(vars(CellMetrics()))

            for _cellwise_records, curr_tag in iter_tag_groups_from_tsv(
                tsv_iterator=tsv_reader
            ):
                # create a metrics aggregator for the current cell
                metric_aggregator = CellMetrics()
                cell_tag = self.compute_metrics(_cellwise_records, metric_aggregator)
                # write a record for each cell
                metric_aggregator.finalize(mitochondrial_genes=self._mitochondrial_gene_ids)
                cell_metrics_output.write(cell_tag, vars(metric_aggregator))


class GatherGeneMetrics(MetricGatherer):

    extra_docs = """
    Notes
    -----
    ``bam_file`` must be sorted by molecule (``UB``), cell (``CB``), and gene (``GE``), where
    molecule varies fastest.

    Examples
    --------
    >>> from sctools.metrics.gatherer import GatherCellMetrics
    >>> import os, tempfile

    >>> # example data
    >>> bam_file = os.path.abspath(__file__) + '../test/data/test.bam'
    >>> temp_dir = tempfile.mkdtemp()
    >>> g = GatherCellMetrics(bam_file=bam_file, output_stem=temp_dir + 'test', compress=True)
    >>> g.extract_metrics()

    See Also
    --------
    GatherGeneMetrics

    """

    __doc__ += extra_docs

    def extract_metrics(self, mode: str = "rb") -> None:
        """Extract gene metrics from self.bam_file

        Parameters
        ----------
        mode : str, optional
            Open mode for self.bam. 'r' -> sam, 'rb' -> bam (default = 'rb').

        """
        # open the files
        with pysam.AlignmentFile(self.bam_file, mode=mode) as bam_iterator, closing(
            MetricCSVWriter(self._output_stem, self._compress)
        ) as gene_metrics_output:

            # write the header
            gene_metrics_output.write_header(vars(GeneMetrics()))

            # break up the bam file into sub-iterators over gene ids
            for gene_iterator, gene_tag in iter_genes(bam_iterator=bam_iterator):
                metric_aggregator = GeneMetrics()

                # in case of multi-genes ignore as in the counting stage
                if gene_tag and len(gene_tag.split(",")) > 1:
                    continue

                # break up gene ids by cell barcodes
                for cell_iterator, cell_tag in iter_cell_barcodes(
                    bam_iterator=gene_iterator
                ):

                    # break up cell barcodes by molecular barcodes
                    for molecule_iterator, molecule_tag in iter_molecule_barcodes(
                        bam_iterator=cell_iterator
                    ):

                        print(gene_tag, cell_tag, molecule_tag)
                        # process the data
                        metric_aggregator.parse_molecule(
                            tags=(gene_tag, cell_tag, molecule_tag),
                            records=molecule_iterator,
                        )

                # write a record for each gene id
                metric_aggregator.finalize()
                gene_metrics_output.write(gene_tag, vars(metric_aggregator))


class GatherGeneMetricsFast(MetricGatherer):

    extra_docs = """
    Notes
    -----
    ``bam_file`` must be sorted by molecule (``UB``), cell (``CB``), and gene (``GE``), where
    molecule varies fastest.

    Examples
    --------
    >>> from sctools.metrics.gatherer import GatherCellMetrics
    >>> import os, tempfile

    >>> # example data
    >>> bam_file = os.path.abspath(__file__) + '../test/data/test.bam'
    >>> temp_dir = tempfile.mkdtemp()
    >>> g = GatherCellMetrics(bam_file=bam_file, output_stem=temp_dir + 'test', compress=True)
    >>> g.extract_metrics()

    See Also
    --------
    GatherGeneMetrics

    """

    __doc__ += extra_docs

    def extract_metrics(self, mode: str = "r") -> None:
        """Extract gene metrics from self.bam_file

        Parameters
        ----------
        mode : str, optional
            Open mode for self.bam. 'r' -> sam, 'rb' -> bam (default = 'rb').

        """

        # open the files note self.bam_file is not always a bam file
        with gzip.open(self.bam_file, mode=mode) if self.bam_file.endswisth('.gz') else \
           open(self.bam_file, mode=mode)  as tsv_reader, closing(
            MetricCSVWriter(self._output_stem, self._compress)
        ) as gene_metrics_output:
            # write the header
            gene_metrics_output.write_header(vars(GeneMetrics()))

            for _cellwise_records, curr_tag in iter_tag_groups_from_tsv(
                tsv_iterator=tsv_reader
            ):
                if curr_tag and len(curr_tag.split(",")) > 1:
                    continue

                metric_aggregator = GeneMetrics()
                gene_tag = self.compute_metrics(_cellwise_records, metric_aggregator)
                # write a record for each gene id
                metric_aggregator.finalize()
                gene_metrics_output.write(gene_tag, vars(metric_aggregator))
