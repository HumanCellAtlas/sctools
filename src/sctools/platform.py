from . import fastq, bam
import argparse


class GenericPlatform:
    """Abstract base class for platform specific metadata"""

    @classmethod
    def attach_barcodes(cls, args=None):
        """command line entrypoint for attaching barcodes to a bamfile."""
        raise NotImplementedError

    @classmethod
    def split_bam(cls, args=None):
        """command line entrypoint for splitting a bamfile into subfiles of equal size.

        prints filenames of chunks to stdout

        :param list args: optional arguments list, for testing (see test/test_entrypoints.py for
          example).
        :return int:
        """
        parser = argparse.ArgumentParser()
        parser.add_argument(
            '-b', '--bamfile', required=True,
            help='input bamfile')
        parser.add_argument(
            '-p', '--output-prefix', required=True,
            help='prefix for output chunks')
        parser.add_argument(
            '-s', '--subfile-size', required=False, default=1000, type=float,
            help='approximate size target for each subfile (in MB)')
        parser.add_argument('-t', '--tag', required=True, help='tag to split bamfile over')
        parser.set_defaults(raise_missing=True)
        parser.add_argument('--drop-missing', action='store_false',
                            help='drop records without tag specified by -t/--tag (default '
                                 'behavior is to raise an exception')
        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()

        filenames = bam.split(
            args.bamfile, args.output_prefix, args.tag, args.subfile_size, args.raise_missing)

        print(' '.join(filenames))
        return 0


class TenXV2(GenericPlatform):

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
    def make_tag_generators(cls, r1, i1=None, whitelist=None):
        """Create tag generators from fastq files.

        Tag generators are iterators that run over fastq records, they extract and yield all of the
        barcodes embedded in each fastq record. For 10x, this means extracting the cell, umi, and
        sample barcode.

        :param str r1: fastq file, read 1
        :param str i1: (optional) fastq file, index read
        :param str whitelist: file containing barcodes -- one barcode perline.
        :return [Generator]: a list of embedded barcodes generators associated with each fastq
          record
        """
        tag_generators = []

        # generator for cell and molecule barcodes
        if whitelist is not None:
            tag_generators.append(fastq.BarcodeGeneratorWithCorrectedCellBarcodes(
                fastq_files=r1,
                embedded_cell_barcode=cls.cell_barcode,
                whitelist=whitelist,
                other_embedded_barcodes=[cls.molecule_barcode],
            ))
        else:
            tag_generators.append(fastq.EmbeddedBarcodeGenerator(
                fastq_files=r1, embedded_barcodes=[cls.cell_barcode, cls.molecule_barcode]))

        # generator for sample barcodes
        if i1 is not None:
            tag_generators.append(fastq.EmbeddedBarcodeGenerator(
                fastq_files=i1, embedded_barcodes=[cls.sample_barcode]))
        return tag_generators

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
            help='read 1 fastq file for a 10x genomics v2 experiment')
        parser.add_argument(
            '--u2', required=True,
            help='unaligned bam containing cDNA fragments. Can be converted from fastq read 2'
                 'using picard FastqToSam')
        parser.add_argument(
            '--i1', default=None,
            help='(optional) i7 index fastq file for a 10x genomics experiment')
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

    @classmethod
    def calculate_cell_metrics(cls, args=None):
        raise NotImplementedError

    @classmethod
    def calculate_gene_metrics(cls, args=None):
        raise NotImplementedError

    @classmethod
    def merge_cell_metrics(cls, args=None):
        raise NotImplementedError

    @classmethod
    def merge_gene_metrics(cls, args=None):
        raise NotImplementedError
