"""
Metric Writers
==============

..currentmodule:: sctools.metrics

This module defines a class to write metrics to csv as the data is generated, cell by cell or gene
by gene. This strategy keeps memory usage low, as no more than a single molecule's worth of sam
records and one cell or gene's worth of metric data are in-memory at a time.

Classes
-------
MetricCSVWriter                 Class to write metrics to file

See Also
--------
sctools.metrics.gatherer
sctools.metrics.aggregator
sctools.metrics.merge

"""
from typing import TextIO, List, Mapping, Any
from numbers import Number
import gzip


class MetricCSVWriter:
    """Writes metric information iteratively to (optionally compressed) csv.

    Parameters
    ----------
    output_stem : str
        File stem for the output file.
    compress : bool, optional
        Whether or not to compress the output file (default = True).

    Methods
    -------
    write_header
        Write the metric header to file.
    write
        Write an array of cell or gene metrics to file.
    close
        Close the metric file.

    """

    def __init__(self, output_stem: str, compress=True):

        # check and fix extension:
        if compress:
            if not output_stem.endswith('.csv.gz'):
                output_stem += '.csv.gz'
        else:
            if not output_stem.endswith('.csv'):
                output_stem += '.csv'
        self._filename: str = output_stem

        # open the file
        if compress:
            self._open_fid: TextIO = gzip.open(self._filename, 'wt')
        else:
            self._open_fid: TextIO = open(self._filename, 'w')
        self._header: List[str] = None

    @property
    def filename(self) -> str:
        """filename with correct suffix added"""
        return self._filename

    def write_header(self, record: Mapping[str, Any]) -> None:
        """Write the metric keys to file, producing the header line of the csv file.

        Parameters
        ----------
        record : Mapping[str, Any]
            Output of ``vars()`` called on an sctools.metrics.aggregator.MetricAggregator instance,
            producing a dictionary of keys to metric values.

        """
        self._header = list(key for key in record.keys() if not key.startswith('_'))
        self._open_fid.write(',' + ','.join(self._header) + '\n')

    def write(self, index: str, record: Mapping[str, Number]) -> None:
        """Write the array of metric values for a cell or gene to file.

        Parameters
        ----------
        index : str
            The name of the cell or gene that these metrics summarize
        record : Mapping[str, Number]
            Output of ``vars()`` called on an sctools.metrics.aggregator.MetricAggregator instance,
            producing a dictionary of keys to metric values.

        """
        ordered_fields = [str(record[k]) for k in self._header]

        # genes and cells can be None, call repr to convert to string when this induces a TypeError
        try:
            self._open_fid.write(index + ',' + ','.join(ordered_fields) + '\n')
        except TypeError:
            index = repr(index)
            self._open_fid.write(index + ',' + ','.join(ordered_fields) + '\n')

    def close(self) -> None:
        """Close the metrics file."""
        self._open_fid.close()
