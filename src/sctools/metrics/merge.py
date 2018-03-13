from typing import Iterable, List
import pandas as pd
import numpy as np


class MergeMetrics:

    def __init__(self, metric_files: Iterable[str], output_file: str):
        """Merges multiple metrics files

        :param metric_files: metrics files to merge
        :param output_file: name for the merged output
        """
        self._metric_files = metric_files
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
        concatenated_frame.to_csv(self._output_file)


class MergeGeneMetrics(MergeMetrics):

    # todo this class needs to do a weighted mean; this is something numpy can do.
    def execute(self) -> None:
        """
        Currently this class returns an incorrect result
        """
        metric_dataframes: List[pd.DataFrame] = [
            pd.read_csv(f, index_col=0) for f in self._metric_files
        ]
        index: pd.Index = metric_dataframes[0].index
        columns: pd.Index = metric_dataframes[0].columns
        averaged_data: np.ndarray = np.mean(metric_dataframes, axis=0)

        # write the data
        pd.DataFrame(data=averaged_data, index=index, columns=columns).to_csv(self._output_file)
