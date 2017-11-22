from . import fastq, bam
import argparse


class PlatformBase:
    """Abstract base class for platform specific metadata"""

    @classmethod
    def attach_barcodes(cls, args=None):
        """command line entrypoint for attaching barcodes to a bamfile."""
        raise NotImplementedError


class TenXV2(PlatformBase):

    # 10x contains three barcodes embedded within sequencing reads. The below objects define the
    # start and end points of those barcodes relative to the start of the sequence, and the
    # GA4GH standard tags that the extracted barcodes should be labeled with in the BAM file.
    cell_barcode = fastq.EmbeddedBarcode(start=0, end=16, quality_tag='CY', sequence_tag='CR')
    molecule_barcode = fastq.EmbeddedBarcode(start=16, end=26, quality_tag='UY', sequence_tag='UR')
    sample_barcode = fastq.EmbeddedBarcode(start=0, end=8, quality_tag='SY', sequence_tag='SR')

    @classmethod
    def tag_bamfile(cls, input_bamfile_name, output_bamfile_name, tag_generators):
        """Extracts tags from fastq files_with_tags, attaches them to input_bamfile_name and writes
        the result to output_bamfile_name

        :param str input_bamfile_name: input bam
        :param str output_bamfile_name: output bam
        :param list tag_generators: list of generators of tags
        """
        bam_tagger = bam.Tagger(input_bamfile_name)
        bam_tagger.tag(output_bamfile_name, tag_generators)

    @classmethod
    def make_tag_generators(cls, r1, i1, whitelist=None):
        """Create tag generators from fastq files.

        Tag generators are iterators that run over fastq records, they extract and yield all of the
        barcodes embedded in each fastq record. For 10x, this means extracting the cell, umi, and
        sample barcode.

        :param str r1:
        :param str i1:
        :param str whitelist:
        :return [Generator]: a list of embedded barcodes generators associated with each fastq
          record
        """
        if whitelist is not None:
            return [
                fastq.BarcodeGeneratorWithCorrectedCellBarcodes(
                    fastq_files=r1,
                    embedded_cell_barcode=cls.cell_barcode,
                    whitelist=whitelist,
                    other_embedded_barcodes=[cls.molecule_barcode],
                ),
                fastq.EmbeddedBarcodeGenerator(
                    fastq_files=i1, embedded_barcodes=[cls.sample_barcode])
            ]
        else:
            return [
                fastq.EmbeddedBarcodeGenerator(
                    fastq_files=r1, embedded_barcodes=[cls.cell_barcode, cls.molecule_barcode]),
                fastq.EmbeddedBarcodeGenerator(
                    fastq_files=i1, embedded_barcodes=[cls.sample_barcode])
            ]

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
        parser.add_argument('-w', '--whitelist', default=None,
                            help='optional cell barcode whitelist. If provided, corrected barcodes '
                                 'will also be output when barcodes are observed within 1ED of a '
                                 'whitelisted barcode')
        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()

        tag_generators = cls.make_tag_generators(args.r1, args.i1, args.whitelist)
        cls.tag_bamfile(args.u2, args.output_bamfile, tag_generators)

        return 0
