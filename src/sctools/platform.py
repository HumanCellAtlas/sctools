"""
Command Line Interface for SC Tools:
====================================

.. currentmodule:: sctools

This module defines the command line interface for SC Tools. Tools are separated into those that
are specific to particular chemistries (e.g. Smart-seq 2) or experimental platforms (e.g. 10x
Genomics v2) and those that are general across any sequencing experiment.

Currently, only general modules and those used for 10x v2 are implemented

Classes
-------
GenericPlatform         Class containing all general command line utilities
TenXV2                  Class containing 10x v2 specific command line utilities

"""

import argparse
from typing import Iterable, List

from sctools import fastq, bam, metrics


class GenericPlatform:
    """Platform-agnostic command line functions available in SC Tools.

    Platform-Agnostic Methods
    -------------------------
    split_bam()
        split a bam file into subfiles of equal size
    calculate_gene_metrics()
        calculate information about genes captured by a sequencing experiment
    calculate_cell_metrics()
        calculate information about cells captured by a sequencing experiment
    merge_gene_metrics()
        merge multiple gene metrics files into a single output
    merge_cell_metrics()
        merge multiple cell metrics files into a single output

    """

    @classmethod
    def split_bam(cls, args: Iterable=None) -> int:
        """Command line entrypoint for splitting a bamfile into subfiles of equal size.

        prints filenames of chunks to stdout

        Parameters
        ----------
        args : Iterable[str], optional
            arguments list, for testing (see test/test_entrypoints.py for example). The default
            value of None, when passed to `parser.parse_args` causes the parser to
            read `sys.argv`

        Returns
        -------
        return_call : 0
            return call if the program completes successfully

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
        parser.add_argument('-t', '--tags', nargs='+',
                            help='tag(s) to split bamfile over. Tags are checked sequentially, '
                                 'and tags after the first are only checked if the first tag is '
                                 'not present.')
        parser.set_defaults(raise_missing=True)
        parser.add_argument('--drop-missing', action='store_false',
                            help='drop records without tag specified by -t/--tag (default '
                                 'behavior is to raise an exception')
        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()

        filenames = bam.split(
            args.bamfile, args.output_prefix, *args.tags,
            approx_mb_per_split=args.subfile_size, raise_missing=args.raise_missing)

        print(' '.join(filenames))
        return 0

    @classmethod
    def calculate_gene_metrics(cls, args: Iterable[str]=None) -> int:
        """Command line entrypoint for calculating gene metrics from a sorted bamfile.

        Writes metrics to .csv

        Parameters
        ----------
        args : Iterable[str], optional
            arguments list, for testing (see test/test_entrypoints.py for example). The default
            value of None, when passed to `parser.parse_args` causes the parser to
            read `sys.argv`

        Returns
        -------
        return_call : 0
            return call if the program completes successfully

        """
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input-bam', required=True, help='Input bam file name.')
        parser.add_argument('-o', '--output-filestem', required=True, help='Output file stem.')

        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()

        gene_metric_gatherer = metrics.gatherer.GatherGeneMetrics(
            args.input_bam, args.output_filestem)
        gene_metric_gatherer.extract_metrics()
        return 0

    @classmethod
    def calculate_cell_metrics(cls, args: Iterable[str]=None) -> int:
        """Command line entrypoint for calculating cell metrics from a sorted bamfile.

        Writes metrics to .csv

        Parameters
        ----------
        args : Iterable[str], optional
            Arguments list, for testing (see test/test_entrypoints.py for example). The default
            value of None, when passed to `parser.parse_args` causes the parser to
            read `sys.argv`

        Returns
        -------
        return_call : 0
            return call if the program completes successfully

        """
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input-bam', required=True, help='Input bam file name.')
        parser.add_argument('-o', '--output-filestem', required=True, help='Output file stem.')

        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()
        cell_metric_gatherer = metrics.gatherer.GatherCellMetrics(
            args.input_bam, args.output_filestem)
        cell_metric_gatherer.extract_metrics()
        return 0

    @classmethod
    def merge_gene_metrics(cls, args: Iterable[str]=None) -> int:
        """Command line entrypoint for merging multiple gene metrics files.

        Merges multiple metrics inputs into a single metrics file that matches the shape and
        order of the generated count matrix.

        Parameters
        ----------
        args : Iterable[str], optional
            Arguments list, for testing (see test/test_entrypoints.py for example). The default
            value of None, when passed to `parser.parse_args` causes the parser to
            read `sys.argv`

        Returns
        -------
        return_call : 0
            return call if the program completes successfully

        """
        parser = argparse.ArgumentParser()
        parser.add_argument('metric_files', nargs='+', help='Input metric files')
        parser.add_argument('-o', '--output-filestem', required=True, help='Output file stem.')

        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()
        merge = metrics.merge.MergeGeneMetrics(args.metric_files, args.output_filestem)
        merge.execute()
        return 0

    @classmethod
    def merge_cell_metrics(cls, args: Iterable[str]=None) -> int:
        """Command line entrypoint for merging multiple cell metrics files.

        Merges multiple metrics inputs into a single metrics file that matches the shape and
        order of the generated count matrix.

        Parameters
        ----------
        args : Iterable[str], optional
            Arguments list, for testing (see test/test_entrypoints.py for example). The default
            value of None, when passed to `parser.parse_args` causes the parser to
            read `sys.argv`

        Returns
        -------
        return_call : 0
            return call if the program completes successfully

        """
        parser = argparse.ArgumentParser()
        parser.add_argument('metric_files', nargs='+', help='Input metric files')
        parser.add_argument('-o', '--output-filestem', required=True, help='Output file stem.')

        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()
        merge = metrics.merge.MergeCellMetrics(args.metric_files, args.output_filestem)
        merge.execute()
        return 0


class TenXV2(GenericPlatform):
    """Command Line Interface for 10x Genomics v2 RNA-sequencing programs

    This class defines several methods that are created as CLI tools when sctools is installed
    (see setup.py)

    Attributes
    ----------
    cell_barcode : fastq.EmbeddedBarcode
        A data class that defines the start and end position of the cell barcode and the tags to
        assign the sequence and quality of the cell barcode
    molecule_barcode : fastq.EmbeddedBarcode
        A data class that defines the start and end position of the molecule barcode and the tags
        to assign the sequence and quality of the molecule barcode
    sample_barcode : fastq.EmbeddedBarcode
        A data class that defines the start and end position of the sample barcode and the tags
        to assign the sequence and quality of the sample barcode

    Methods
    -------
    attach_barcodes()
        Attach barcodes from the forward (r1) and optionally index (i1) fastq files to the reverse
        (r2) bam file

    """
    # 10x contains three barcodes embedded within sequencing reads. The below objects define the
    # start and end points of those barcodes relative to the start of the sequence, and the
    # GA4GH standard tags that the extracted barcodes should be labeled with in the BAM file.
    cell_barcode = fastq.EmbeddedBarcode(start=0, end=16, quality_tag='CY', sequence_tag='CR')
    molecule_barcode = fastq.EmbeddedBarcode(start=16, end=26, quality_tag='UY', sequence_tag='UR')
    sample_barcode = fastq.EmbeddedBarcode(start=0, end=8, quality_tag='SY', sequence_tag='SR')

    @classmethod
    def _tag_bamfile(
            cls,
            input_bamfile_name: str,
            output_bamfile_name: str,
            tag_generators: Iterable[fastq.EmbeddedBarcodeGenerator]) -> None:
        """Adds tags from fastq file(s) to a bam file.

        Attaches tags extracted from fastq files by `tag_generators`, attaches them to records from
        `input_bamfile_name`, and writes the result to `output_bamfile_name`

        Parameters
        ----------
        input_bamfile_name : str
            input bam
        output_bamfile_name : str
            output bam
        tag_generators : Iterable[fastq.EmbeddedBarcodeGenerator]
            Iterable of generators that yield barcodes from fastq files

        """
        bam_tagger = bam.Tagger(input_bamfile_name)
        bam_tagger.tag(output_bamfile_name, tag_generators)

    @classmethod
    def _make_tag_generators(
            cls, r1, i1=None, whitelist=None) -> List[fastq.EmbeddedBarcodeGenerator]:
        """Create tag generators from fastq files.

        Tag generators are iterators that run over fastq records, they extract and yield all of the
        barcodes embedded in each fastq record. For 10x, this means extracting the cell, umi, and
        optionally, the sample barcode.

        Parameters
        ----------
        r1 : str
            forward fastq file
        i1 : str, optional
            index fastq file
        whitelist : str, optional
            A file that contains a list of acceptable cell barcodes

        Returns
        -------
        tag_generators, List[EmbeddedBarcodeGenerator]
            EmbeddedBarcodeGenerators containing barcodes from 10x fastq records

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
        """Command line entrypoint for attaching barcodes to a bamfile.

        Parameters
        ----------
        args : Iterable[str], optional
            arguments list, for testing (see test/test_entrypoints.py for example). The default
            value of None, when passed to `parser.parse_args` causes the parser to
            read `sys.argv`

        Returns
        -------
        return_call : 0
            return call if the program completes successfully

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
        tag_generators = cls._make_tag_generators(args.r1, args.i1, args.whitelist)
        cls._tag_bamfile(args.u2, args.output_bamfile, tag_generators)

        return 0
