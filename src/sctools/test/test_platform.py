import os
import tempfile
import pysam

from .. import platform

data_dir = os.path.split(__file__)[0] + "/data/"


def test_attach_barcodes():
    """High-level test of the AttachBarcodes command"""

    temp_dir_name = tempfile.mkdtemp()

    # Construct cli arguments to pass to the command
    temp_output_bam = temp_dir_name + "output.bam"

    args = [
        "--r1",
        data_dir + "test_r1.fastq",
        "--u2",
        data_dir + "test_r2.bam",
        "--i1",
        data_dir + "test_i1.fastq",
        "--o",
        temp_output_bam,
        "--sample-barcode-start-pos",
        "0",
        "--sample-barcode-length",
        "8",
        "--cell-barcode-start-pos",
        "0",
        "--cell-barcode-length",
        "16",
        "--molecule-barcode-start-pos",
        "16",
        "--molecule-barcode-length",
        "4",
    ]

    platform.BarcodePlatform.attach_barcodes(args)

    with pysam.AlignmentFile(temp_output_bam, "rb", check_sq=False) as samfile:
        for read in samfile:
            tag_cr = read.get_tag("CR")
            tag_cy = read.get_tag("CY")
            tag_ur = read.get_tag("UR")
            tag_uy = read.get_tag("UY")
            tag_sr = read.get_tag("SR")
            tag_sy = read.get_tag("SY")
            assert len(tag_cr) == 16
            assert len(tag_cy) == 16
            assert len(tag_ur) == 4
            assert len(tag_uy) == 4
            assert len(tag_sr) == 8
            assert len(tag_sy) == 8
