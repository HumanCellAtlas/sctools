from . import fastq, bam
import argparse


class PlatformBase:
    """Abstract base class for platform specific metadata"""

    @classmethod
    def get_tags(cls, sequencing_read):
        raise NotImplementedError

    @classmethod
    def tag_bamfile(cls, input_bamfile_name, output_bamfile_name, files_with_tags):
        """Extracts tags from fastq files_with_tags, attaches them to input_bamfile_name and writes
        the result to output_bamfile_name

        :param str input_bamfile_name: input bam
        :param str output_bamfile_name: output bam
        :param dict files_with_tags: dict mapping {tags: filenames}
        """
        tag_generators = []
        for k, v in files_with_tags.items():
            tag_generators.append(fastq.EmbeddedBarcodeGenerator(cls.get_tags(k), files=v))

        bam_tagger = bam.Tagger(input_bamfile_name)
        bam_tagger.tag(output_bamfile_name, tag_generators)

    @classmethod
    def attach_barcodes(cls, args=None):
        """command line entrypoint for attaching barcodes to a bamfile."""
        raise NotImplementedError


class TenXV2(PlatformBase):

    cell_barcode = fastq.EmbeddedBarcode(start=0, end=16, quality_tag='CY', sequence_tag='CR')
    molecule_barcode = fastq.EmbeddedBarcode(start=16, end=24, quality_tag='UY', sequence_tag='UR')
    sample_barcode = fastq.EmbeddedBarcode(start=0, end=8, quality_tag='SY', sequence_tag='SR')

    _tags = {
        'r1': (cell_barcode, molecule_barcode),
        'i1': (sample_barcode,)
    }

    @classmethod
    def get_tags(cls, sequencing_read):
        return cls._tags[sequencing_read]

    @classmethod
    def attach_barcodes(cls, args=None):
        """command line entrypoint for attaching barcodes to a bamfile.
        :param list args: optional arguments list, for testing (see test/test_entrypoints.py for
          example).
        :return int:
        """

        parser = argparse.ArgumentParser()
        parser.add_argument(
            '--r1', required=True,
            help='read 1 fastq record for a 10x genomics experiment')
        parser.add_argument(
            '--i1', required=True, help='i7 fastq record for a 10x genomics experiment')
        parser.add_argument(
            '--u2', required=True,
            help='unaligned read-2 bam containing genomic information. Can be converted'
                 'using picard FastqToSam')
        parser.add_argument('-o', '--output-bamfile', required=True,
                            help='filename for tagged bam')
        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()

        cls.tag_bamfile(args.u2, args.output_bamfile, {'r1': args.r1, 'i1': args.i1})

        return 0
