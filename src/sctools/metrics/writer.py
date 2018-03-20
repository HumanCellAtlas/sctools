from typing import TextIO, List, Mapping, Any
from numbers import Number


class MetricCSVWriter:

    def __init__(self, output_stem: str):
        self._filename: str = output_stem if output_stem.endswith('.csv') else output_stem + '.csv'
        self._open_fid: TextIO = open(self._filename, 'w')
        self._header: List[str] = None

    def filename(self) -> str:
        return self._filename

    def write_header(self, record: Mapping[str, Any]) -> None:
        self._header = list(key for key in record.keys() if not key.startswith('_'))
        self._open_fid.write(',' + ','.join(self._header) + '\n')

    def write(self, index: str, record: Mapping[str, Number]) -> None:
        ordered_fields = [str(record[k]) for k in self._header]
        self._open_fid.write(index + ',' + ','.join(ordered_fields) + '\n')

    def close(self) -> None:
        self._open_fid.close()
