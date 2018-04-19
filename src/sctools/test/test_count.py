"""
Testing for Count Matrix Construction
=====================================

        bam_file: str,
        annotation_file: str,
        cell_barcode_tag: str='CB',
        molecule_barcode_tag: str='UB',
        gene_id_tag: str='GE',
        open_mode: str='rb',

"""
import os
import tempfile
from sctools.count import bam_to_count

# set the input and output directories, using a tempdir to automatically clean up generated files
_data_dir = os.path.split(__file__)[0] + '/data'
_test_dir = tempfile.mkdtemp()
os.makedirs(_test_dir, exist_ok=True)

_bam_file = _data_dir + '/test.bam'
_gtf_file = _data_dir + '/test.gtf.gz'


def test_construct_count_matrix():
    count_matrix = bam_to_count(_bam_file, _gtf_file)
    assert count_matrix.shape == (0, 16)

# todo need to add test data that includes counts
