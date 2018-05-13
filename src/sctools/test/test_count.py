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

import numpy as np
import scipy.sparse as sp

from sctools.count import CountMatrix

# set the input and output directories, using a tempdir to automatically clean up generated files
_data_dir = os.path.split(__file__)[0] + '/data'
_test_dir = tempfile.mkdtemp()
os.makedirs(_test_dir, exist_ok=True)

_bam_file = _data_dir + '/cell-sorted2.bam'
_gtf_file = _data_dir + '/chr1.gtf.gz'


# todo this test should be made faster
def test_count_matrix():
    count_matrix = CountMatrix.from_bam(_bam_file, _gtf_file)
    assert count_matrix._matrix.shape == (734, 3132)

    # check save + load
    count_matrix.save(_test_dir + '/test_count_matrix')
    count_matrix_2 = CountMatrix.load(_test_dir + '/test_count_matrix')

    # check matrices are the same
    csr1: sp.csr_matrix = count_matrix._matrix
    csr2: sp.csr_matrix = count_matrix_2._matrix
    assert np.allclose(csr1.indices, csr2.indices)
    assert np.allclose(csr1.indptr, csr2.indptr)
    assert np.allclose(csr1.data, csr2.data)
