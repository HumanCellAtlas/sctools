from typing import List, Sequence
import pandas as pd
import numpy as np


class MergeMetrics:

    def __init__(self, metric_files: Sequence[str], output_file: str):
        """Merges multiple metrics files

        :param metric_files: metrics files to merge
        :param output_file: name for the merged output
        """
        self._metric_files = metric_files
        if not output_file.endswith('.csv.gz'):
            output_file += '.csv.gz'
        self._output_file = output_file

    def execute(self) -> None:
        raise NotImplementedError  # merge the metrics


class MergeCellMetrics(MergeMetrics):

    # todo test me
    def execute(self) -> None:
        """
        shutil would be more efficient, but we need a way to get rid of the header line. Bash
        might be faster and more efficient than the python solution
        """
        metric_dataframes: List[pd.DataFrame] = [
            pd.read_csv(f, index_col=0) for f in self._metric_files
        ]
        concatenated_frame: pd.DataFrame = pd.concat(metric_dataframes, axis=0)
        concatenated_frame.to_csv(self._output_file, compression='gzip')


class MergeGeneMetrics(MergeMetrics):

    # todo this class needs to do a weighted mean; this is something numpy can do.
    def execute(self) -> None:
        """
        Genes are expected to be found in multiple shards. Thus, the metric values must be
        summed or averaged over, depending on how they were originally calculated.

        This function exectures the merge and writes the resulting data to self.output_file
        """

        # different columns in the gene dataset need to be merged in different ways. Count data
        # must be summed, while averages need to have a (weighted) average taken

        columns_to_sum = [
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

        sum_operations = {c: 'sum' for c in columns_to_sum}

        def weighted_average(data_frame):
            """
            :param pd.DataFrame data_frame: input dataframe to reduce with a weighted average
            :return pd.Series: reduced result
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
                {c: np.average(data_frame[c], weights=weights) for c in columns_to_average_by_read})

        def recalculate_operation(x):
            return pd.DataFrame(
                data={
                    'reads_per_molecule': x['n_reads'] / x['n_molecules'],
                    'fragments_per_molecule': x['n_fragments'] / x['n_molecules'],
                    'reads_per_fragment': x['n_reads'] / x['n_fragments']
                }
            )

        # now merge each subsequent dataframe into the first one
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
