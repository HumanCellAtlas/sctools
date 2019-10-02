import argparse
import os
import tempfile
import pysam

from .. import platform

data_dir = os.path.split(__file__)[0] + '/data/'

def test_attach_barcodes():
    """High-level test of the AttachBarcodes command"""

    temp_dir_name = tempfile.mkdtemp()

    # Construct cli arguments to pass to the command
    temp_output_bam = temp_dir_name + 'output.bam'
    print(temp_output_bam);
    args = [
        "--r1", data_dir + 'test_r1.fastq',
        "--u2", data_dir + 'test_r2.bam',
        "--i1", data_dir + 'test_i1.fastq',
        "--o", temp_output_bam,
        "--sample-barcode-start-pos", "0",
        "--sample-barcode-length", "8",
        "--cell-barcode-start-pos", "0",
        "--cell-barcode-length", "16",
        "--molecule-barcode-start-pos", "16",
        "--molecule-barcode-length", "12"
    ]

    platform.BarcodePlatform.attach_barcodes(args)

    with pysam.AlignmentFile(temp_output_bam, "rb", check_sq=False) as samfile:
        for read in samfile:
            tag_CR = read.get_tag('CR')
            tag_CY = read.get_tag('CY')
            tag_UR = read.get_tag('UR')
            tag_UY = read.get_tag('UY')
            tag_SR = read.get_tag('SR')
            tag_SY = read.get_tag('SY')
            assert len(tag_CR) == 16
            assert len(tag_CY) == 16
            assert len(tag_UR) == 11
            assert len(tag_UY) == 11
            assert len(tag_SR) == 8
            assert len(tag_SY) == 8
