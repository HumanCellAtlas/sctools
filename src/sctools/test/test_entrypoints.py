from .. import platform
import os
import pysam

data_dir = os.path.split(__file__)[0] + '/data/'


def test_Attach10XBarcodes_entrypoint():
    args = [
        '--r1', data_dir + 'test_r1.fastq',
        '--i1', data_dir + 'test_i7.fastq',
        '--u2', data_dir + 'test_r2.bam',
        '--output-bamfile', 'test_tagged_bam.bam']

    rc = platform.TenXV2.attach_barcodes(args)
    assert rc == 0
    with pysam.AlignmentFile('test_tagged_bam.bam', 'rb', check_sq=False) as f:
        for alignment in f:
            # each alignment should now have a tag, and that tag should be a string
            assert isinstance(alignment.get_tag('CY'), str)
            assert isinstance(alignment.get_tag('CR'), str)
            assert isinstance(alignment.get_tag('UY'), str)
            assert isinstance(alignment.get_tag('UR'), str)
            assert isinstance(alignment.get_tag('SY'), str)
            assert isinstance(alignment.get_tag('SR'), str)
    os.remove('test_tagged_bam.bam')  # clean up
