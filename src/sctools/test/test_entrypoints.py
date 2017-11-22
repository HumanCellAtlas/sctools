import os
import glob
import pysam
from .. import platform

data_dir = os.path.split(__file__)[0] + '/data/'


def test_Attach10XBarcodes_entrypoint():
    args = [
        '--r1', data_dir + 'test_r1.fastq',
        '--i1', data_dir + 'test_i7.fastq',
        '--u2', data_dir + 'test.bam',
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


def test_Attach10XBarcodes_entrypoint_with_whitelist():
    args = [
        '--r1', data_dir + 'test_r1.fastq',
        '--i1', data_dir + 'test_i7.fastq',
        '--u2', data_dir + 'test.bam',
        '--output-bamfile', 'test_tagged_bam.bam',
        '--whitelist', data_dir + '1k-august-2016.txt']

    return_call = platform.TenXV2.attach_barcodes(args)
    assert return_call == 0
    success = False
    with pysam.AlignmentFile('test_tagged_bam.bam', 'rb', check_sq=False) as f:
        for alignment in f:
            if alignment.has_tag('CB'):
                success = True
            # each alignment should now have a tag, and that tag should be a string
            assert isinstance(alignment.get_tag('CY'), str)
            assert isinstance(alignment.get_tag('CR'), str)
            assert isinstance(alignment.get_tag('UY'), str)
            assert isinstance(alignment.get_tag('UR'), str)
            assert isinstance(alignment.get_tag('SY'), str)
            assert isinstance(alignment.get_tag('SR'), str)
    assert success
    os.remove('test_tagged_bam.bam')  # clean up


def test_split_bam():
    tag_args = [
        '--r1', data_dir + 'test_r1.fastq',
        '--i1', data_dir + 'test_i7.fastq',
        '--u2', data_dir + 'test.bam',
        '--output-bamfile', 'test_tagged_bam.bam',
        '--whitelist', data_dir + '1k-august-2016.txt']

    platform.TenXV2.attach_barcodes(tag_args)

    split_args = [
        '--bamfile', 'test_tagged_bam.bam',
        '--output-prefix', 'test_tagged',
        '--subfile-size', '0.005',
        '--tag', 'CB']

    return_call = platform.GenericPlatform.split_bam(split_args)
    assert return_call == 0

    for f in glob.glob('test_tagged*'):
        os.remove(f)

