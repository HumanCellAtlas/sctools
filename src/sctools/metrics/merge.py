"""
Merge Sequence Metrics
======================

..currentmodule:: sctools.metrics

This module defines classes to merge multiple metrics files that have been gathered from bam files
containing disjoint sets of cells. This is a common use pattern, as sequencing datasets are often
chunked to enable horizontal scaling using scatter-gather patterns.

Classes
-------
MergeMetrics                    Merge Metrics base class
MergeCellMetrics                Class to merge cell metrics
MergeGeneMetrics                Class to merge gene metrics

See Also
--------
sctools.metrics.gatherer
sctools.metrics.aggregator
sctools.metrics.writer

"""

from typing import List, Sequence

import pandas as pd
import numpy as np


class MergeMetrics:
    """Merges multiple metrics files into a single gzip compressed csv file

    Parameters
    ----------
    metric_files : Sequence[str]
        metrics files to merge
    output_file : str
        file name for the merged output

    Methods
    -------
    execute
        merge metrics files
        # todo this should probably be wrapped into __init__ to make this more like a function

    """

    def __init__(self, metric_files: Sequence[str], output_file: str):
        self._metric_files = metric_files
        if not output_file.endswith('.csv.gz'):
            output_file += '.csv.gz'
        self._output_file = output_file

    def execute(self) -> None:
        raise NotImplementedError  # merge the metrics


class MergeCellMetrics(MergeMetrics):
    def execute(self) -> None:
        """Concatenate input cell metric files

        Since bam files that metrics are calculated from contain disjoint sets of cells, cell
        metrics can simply be concatenated together.

        """
        metric_dataframes: List[pd.DataFrame] = [
            pd.read_csv(f, index_col=0) for f in self._metric_files
        ]
        concatenated_frame: pd.DataFrame = pd.concat(metric_dataframes, axis=0)
        concatenated_frame.to_csv(self._output_file, compression='gzip')


class MergeGeneMetrics(MergeMetrics):
    def execute(self) -> None:
        """Merge input gene metric files

        The bam files that metrics are calculated from contain disjoint sets of cells, each
        of which can measure the same genes.
        As a result, the metric values must be summed (count based metrics) averaged over
        (fractional, averge, or variance metrics) or recalculated (metrics that depend on other
        metrics).

        """

        count_data_to_sum = [
            'n_reads',
            'noise_reads',
            'perfect_molecule_barcodes',
            'reads_mapped_exonic',
            'reads_mapped_intronic',
            'reads_mapped_utr',
            'reads_mapped_uniquely',
            'reads_mapped_multiple',
            'duplicate_reads',
            'spliced_reads',
            'antisense_reads',
            'n_molecules',
            'n_fragments',
            'fragments_with_single_read_evidence',
            'molecules_with_single_read_evidence',
            'number_cells_detected_multiple',
            'number_cells_expressing',
        ]

        sum_operations = {c: 'sum' for c in count_data_to_sum}

        def weighted_average(data_frame: pd.DataFrame) -> pd.Series:
            """Calculate the average of each metric, weighted by number of reads per chunk

            Parameters
            ----------
            data_frame : pd.DataFrame
              chunks x metrics data frame

            Returns
            -------
            weighted_average_metrics : pd.Series
                The average of each metric across chunks, weighted by the number of reads per chunk

            """
            weights = data_frame['n_reads'].values

            columns_to_average_by_read = [
                'molecule_barcode_fraction_bases_above_30_mean',
                'molecule_barcode_fraction_bases_above_30_variance',
                'genomic_reads_fraction_bases_quality_above_30_mean',
                'genomic_reads_fraction_bases_quality_above_30_variance',
                'genomic_read_quality_mean',
                'genomic_read_quality_variance',
            ]

            return pd.Series(
                {
                    c: np.average(data_frame[c], weights=weights)
                    for c in columns_to_average_by_read
                }
            )

        def recalculate_operation(data_frame) -> pd.DataFrame:
            """Recalculate metrics that are dependent on other metric values

            Other metrics should be merged before this function is executed

            Parameters
            ----------
            data_frame : pd.DataFrame
                chunks x metrics data frame

            Returns
            -------
            recalculated_metrics : pd.DataFrame
                data frame containing recalculated metrics

            """
            return pd.DataFrame(
                data={
                    'reads_per_molecule': data_frame['n_reads']
                    / data_frame['n_molecules'],
                    'fragments_per_molecule': data_frame['n_fragments']
                    / data_frame['n_molecules'],
                    'reads_per_fragment': data_frame['n_reads']
                    / data_frame['n_fragments'],
                }
            )

        # pick one file as a nucleus and merge each subsequent dataframe into it
        nucleus = pd.read_csv(self._metric_files[0], index_col=0)
        for filename in self._metric_files[1:]:
            leaf = pd.read_csv(filename, index_col=0)

            # concatenate this leaf with the nucleus metrics file
            concatenated = pd.concat([nucleus, leaf], axis=0)

            # group all duplicate gene names together
            grouped = concatenated.groupby(level=0, axis=0)

            # execute the merging operations
            summed_columns = grouped.agg(sum_operations)
            averaged_columns = grouped.apply(weighted_average)

            # stitch the columns back together, add the metrics that need to be recalculated
            merged = pd.concat([summed_columns, averaged_columns], axis=1)
            recalculated_columns = recalculate_operation(merged)
            merged = pd.concat([merged, recalculated_columns], axis=1)

            # set as nucleus and continue
            nucleus = merged

        # write the data
        nucleus.to_csv(self._output_file, compression='gzip')
