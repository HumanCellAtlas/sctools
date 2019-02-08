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
from typing import Iterable, List, Dict, Optional, Sequence
from itertools import chain

import pysam
from sctools import fastq, bam, metrics, count, consts, gtf, groups


class GenericPlatform:
    """Platform-agnostic command line functions available in SC Tools.

    Platform-Agnostic Methods
    -------------------------
    tag_sort_bam():
        sort a bam file by zero or more tags and then by queryname
    verify_bam_sort():
        verifies whether bam file is correctly sorted by given list of zero or more tags, then queryname
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
    bam_to_count()
        construct a compressed sparse row count file from a tagged, aligned bam file
    merge_count_matrices()
        merge multiple csr-format count matrices into a single csr matrix
    group_qc_outputs()
        aggregate Picard, HISAT2 and RSME QC statisitics
    """

    @classmethod
    def tag_sort_bam(cls, args: Iterable=None) -> int:
        """Command line entrypoint for sorting a bam file by zero or more tags, followed by queryname.

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
        description = 'Sorts bam by list of zero or more tags, followed by query name'
        parser = argparse.ArgumentParser(description=description)
        parser.add_argument(
            '-i', '--input_bam', required=True,
            help='input bamfile')
        parser.add_argument(
            '-o', '--output_bam', required=True,
            help='output bamfile')
        parser.add_argument('-t', '--tags', nargs='+', action='append',
                            help='tag(s) to sort by, separated by space, e.g. -t CB GE UB')
        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()

        tags = cls.get_tags(args.tags)
        with pysam.AlignmentFile(args.input_bam, 'rb') as f:
            header = f.header
            records = f.fetch(until_eof=True)
            sorted_records = bam.sort_by_tags_and_queryname(records, tags)
        with pysam.AlignmentFile(args.output_bam, 'wb', header=header) as f:
            for record in sorted_records:
                f.write(record)

        return 0

    @classmethod
    def verify_bam_sort(cls, args: Iterable=None) -> int:
        """Command line entrypoint for verifying bam is properly sorted by zero or more tags, followed by queryname.

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
        description = 'Verifies whether bam is sorted by the list of zero or more tags, followed by query name'
        parser = argparse.ArgumentParser(description=description)
        parser.add_argument(
            '-i', '--input_bam', required=True,
            help='input bamfile')
        parser.add_argument('-t', '--tags', nargs='+', action='append',
                            help='tag(s) to use to verify sorting, separated by space, e.g. -t CB GE UB')
        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()

        tags = cls.get_tags(args.tags)
        with pysam.AlignmentFile(args.input_bam, 'rb') as f:
            aligned_segments = f.fetch(until_eof=True)
            sortable_records = (bam.TagSortableRecord.from_aligned_segment(r, tags) for r in aligned_segments)
            bam.verify_sort(sortable_records, tags)

        print('{0} is correctly sorted by {1} and query name'.format(args.input_bam, tags))
        return 0

    @classmethod
    def get_tags(cls, raw_tags: Optional[Sequence[str]]) -> Iterable[str]:
        if raw_tags is None:
            raw_tags = []
        # Flattens into single list when tags specified like -t A -t B -t C
        return [t for t in chain.from_iterable(raw_tags)]

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

    @classmethod
    def bam_to_count_matrix(cls, args: Iterable[str]=None) -> int:
        """Command line entrypoint for constructing a count matrix from a tagged bam file.

        Constructs a count matrix from an aligned bam file sorted by cell barcode, molecule
        barcode, and gene id.

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
        parser.set_defaults(
            cell_barcode_tag=consts.CELL_BARCODE_TAG_KEY,
            molecule_barcode_tag=consts.MOLECULE_BARCODE_TAG_KEY,
            gene_name_tag=consts.GENE_NAME_TAG_KEY
        )
        parser.add_argument('-b', '--bam-file', help='input_bam_file', required=True)
        parser.add_argument(
            '-o', '--output-prefix', help='file stem for count matrix', required=True)
        parser.add_argument(
            '-a', '--gtf-annotation-file', required=True,
            help='gtf annotation file that bam_file was aligned against')
        parser.add_argument(
            '-c', '--cell-barcode-tag',
            help=f'tag that identifies the cell barcode (default = {consts.CELL_BARCODE_TAG_KEY})')
        parser.add_argument(
            '-m', '--molecule-barcode-tag',
            help=f'tag that identifies the molecule barcode (default = {consts.MOLECULE_BARCODE_TAG_KEY})')
        parser.add_argument(
            '-g', '--gene-id-tag',
            help=f'tag that identifies the gene name (default = {consts.GENE_NAME_TAG_KEY})')

        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()

        # assume bam file unless the file explicitly has a sam suffix
        open_mode = 'r' if args.bam_file.endswith('.sam') else 'rb'

        # load gene names from the annotation file
        gene_name_to_index: Dict[str, int] = gtf.extract_gene_names(args.gtf_annotation_file)

        matrix = count.CountMatrix.from_sorted_tagged_bam(
            bam_file=args.bam_file,
            gene_name_to_index=gene_name_to_index,
            cell_barcode_tag=args.cell_barcode_tag,
            molecule_barcode_tag=args.molecule_barcode_tag,
            gene_name_tag=args.gene_name_tag,
            open_mode=open_mode
        )
        matrix.save(args.output_prefix)

        return 0

    @classmethod
    def merge_count_matrices(cls, args: Iterable[str]=None) -> int:
        """Command line entrypoint for constructing a count matrix from a tagged bam file.

        Constructs a count matrix from an aligned bam file sorted by cell barcode, molecule
        barcode, and gene id.

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
        parser.add_argument('-i', '--input-prefixes', nargs='+',
                            help='prefix for count matrices to be concatenated. e.g. test_counts '
                                 'for test_counts.npz, test_counts_col_index.npy, and test_counts_'
                                 'row_index.npy')
        parser.add_argument('-o', '--output-stem', help='file stem for merged csr matrix',
                            required=True)

        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()

        count_matrix = count.CountMatrix.merge_matrices(args.input_prefixes)
        count_matrix.save(args.output_stem)

        return 0

    @classmethod
    def group_qc_outputs(cls, args: Iterable[str]=None) -> int:
        """Commandline entrypoint for parsing picard metrics files, hisat2 and rsem statistics log files.
        Parameters
        ----------
        args:
            file_names: array of files
            output_name: prefix of output file name.
            metrics_type: Picard, PicardTable, HISAT2, RSEM and Core.
        Returns
        ----------
        return: 0 
            return if the program completes successfully.
        """
        parser = argparse.ArgumentParser()
        parser.add_argument(
                "-f",
                "--file_names",
                dest="file_names",
                nargs='+',
                required=True,
                help="a list of files to be parsed out.")
        parser.add_argument(
                "-o",
                "--output_name",
                dest="output_name",
                required=True,
                help="The output file name")
        parser.add_argument(
                "-t",
                "--metrics_type",
                dest="metrics_type",
                choices=['Picard', 'PicardTable', 'Core', 'HISAT2', 'RSEM'],
                required=True,
                help="a list of string to represent metrics types,such Picard, PicardTable, HISAT2,RSEM, Core")

        if args is not None:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()

        if args.metrics_type == "Picard":
            groups.write_aggregated_picard_metrics_by_row(args.file_names, args.output_name)
        elif args.metrics_type == "PicardTable":
            groups.write_aggregated_picard_metrics_by_table(args.file_names, args.output_name)
        elif args.metrics_type == "Core":
            groups.write_aggregated_qc_metrics(args.file_names, args.output_name)
        elif args.metrics_type == "HISAT2":
            groups.parse_hisat2_log(args.file_names, args.output_name)
        elif args.metrics_type == "RSEM":
            groups.parse_rsem_cnt(args.file_names, args.output_name)
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
    cell_barcode = fastq.EmbeddedBarcode(
        start=0,
        end=16,
        quality_tag=consts.QUALITY_CELL_BARCODE_TAG_KEY,
        sequence_tag=consts.RAW_CELL_BARCODE_TAG_KEY)
    molecule_barcode = fastq.EmbeddedBarcode(
        start=16,
        end=26,
        quality_tag=consts.QUALITY_MOLECULE_BARCODE_TAG_KEY,
        sequence_tag=consts.RAW_MOLECULE_BARCODE_TAG_KEY)
    sample_barcode = fastq.EmbeddedBarcode(
        start=0,
        end=8,
        quality_tag=consts.QUALITY_SAMPLE_BARCODE_TAG_KEY,
        sequence_tag=consts.RAW_SAMPLE_BARCODE_TAG_KEY)

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


class Attach(GenericPlatform):

    cell_barcode = None
    molecule_barcode = None
    sample_barcode = None

    @classmethod
    def _get_barcode(cls, barcode_start_pos, barcode_length, barcode_quality_tag, barcode_sequence_tag):
        barcode = {"start": barcode_start_pos,
                   "end": barcode_length,
                   "quality_tag": barcode_quality_tag,
                   "sequence_tag": barcode_sequence_tag}
        return fastq.EmbeddedBarcode(**barcode)

    @classmethod
    def _validate_barcode_input(cls, given_value, min_value):
        if given_value < min_value:
            raise argparse.ArgumentTypeError("Invalid barcode length/position")
        return given_value

    @classmethod
    def _validate_barcode_start_pos(cls, given_value):
        return cls._validate_barcode_input(int(given_value), 0)

    @classmethod
    def _validate_barcode_length(cls, given_value):
        return cls._validate_barcode_input(int(given_value), 1)

    @classmethod
    def _tag_bamfile(cls,
                     input_bamfile_name: str,
                     output_bamfile_name: str,
                     tag_generators: Iterable[fastq.EmbeddedBarcodeGenerator]) -> None:
        bam_tagger = bam.Tagger(input_bamfile_name)
        bam_tagger.tag(output_bamfile_name, tag_generators)

    @classmethod
    def _make_tag_generators(cls, r1, i1=None, whitelist=None) -> List[fastq.EmbeddedBarcodeGenerator]:
        tag_generators = []
        barcode_args = {"fastq_files": r1}

        if i1:
            barcode_args["embedded_barcodes"] = [cls.sample_barcode]
            tag_generators.append(fastq.EmbeddedBarcodeGenerator(**barcode_args))

        if whitelist:
            barcode_args["whitelist"] = whitelist
            if cls.cell_barcode:
                barcode_args["embedded_cell_barcode"] = cls.cell_barcode
            if cls.molecule_barcode:
                barcode_args["other_embedded_barcodes"] = cls.molecule_barcode

        else:
            if cls.cell_barcode and cls.molecule_barcode:
                barcode_args["embedded_barcodes"] = [cls.cell_barcode,
                                                     cls.molecule_barcode]
            elif cls.cell_barcode:
                barcode_args["embedded_barcodes"] = [cls.cell_barcode]
            elif cls.molecule_barcode:
                barcode_args["embedded_barcodes"] = [cls.molecule_barcode]
            tag_generators.append(fastq.EmbeddedBarcodeGenerator(**barcode_args))

        return tag_generators

    @classmethod
    def attach_barcodes(cls, args=None):
        parser = argparse.ArgumentParser()
        parser.add_argument('--r1',
                            required=True,
                            help='read 1 fastq file for a 10x genomics v2 experiment')
        parser.add_argument('--u2',
                            required=True,
                            help='unaligned bam containing cDNA fragments. Can be converted from fastq read 2'
                                 'using picard FastqToSam')
        parser.add_argument('-o',
                            '--output-bamfile',
                            required=True,
                            help='filename for tagged bam')
        parser.add_argument('-w',
                            '--whitelist',
                            default=None,
                            help='optional cell barcode whitelist. If provided, corrected barcodes '
                                 'will also be output when barcodes are observed within 1ED of a '
                                 'whitelisted barcode')
        parser.add_argument('--i1',
                            default=None,
                            help='(optional) i7 index fastq file for a 10x genomics experiment')
        parser.add_argument("--sample-barcode-start-position",
                            dest="sample_barcode_start_pos",
                            default=None,
                            help='the user defined start position (base pairs) of the sample barcode',
                            type=cls._validate_barcode_start_pos)
        parser.add_argument("--sample-barcode-length",
                            dest="sample_barcode_length",
                            default=None,
                            help='the user defined length (base pairs) of the sample barcode',
                            type=cls._validate_barcode_length)
        parser.add_argument("--cell-barcode-start-position",
                            dest="cell_barcode_start_pos",
                            default=None,
                            help='the user defined start position, in base pairs, of the cell barcode',
                            type=cls._validate_barcode_start_pos)
        parser.add_argument("--cell-barcode-length",
                            dest="cell_barcode_length",
                            default=None,
                            help='the user defined length, in base pairs, of the cell barcode',
                            type=cls._validate_barcode_length)
        parser.add_argument("--molecule-barcode-start-position",
                            dest="molecule_barcode_start_pos",
                            default=None,
                            help='the user defined start position, in base pairs, of the molecule barcode '
                                 '(must be not overlap cell barcode if cell barcode is provided)',
                            type=cls._validate_barcode_start_pos)
        parser.add_argument("--molecule-barcode-length",
                            dest="molecule_barcode_length",
                            default=None,
                            help='the user defined length, in base pairs, of the molecule barcode',
                            type=cls._validate_barcode_length)
        if args:
            args = parser.parse_args(args)
        else:
            args = parser.parse_args()

        if ((bool(args.cell_barcode_start_pos) or args.cell_barcode_start_pos == 0)
                != bool(args.cell_barcode_length) or
            (bool(args.molecule_barcode_start_pos) or args.molecule_barcode_start_pos == 0)
                != bool(args.molecule_barcode_length) or
            (bool(args.sample_barcode_start_pos) or args.sample_barcode_start_pos == 0)
                != bool(args.sample_barcode_length)):
            argparse.ArgumentError("Invalid barocde pos/length arguments, barcode start pos and barcode length must be specified together")
        if args.i1 is None and args.sample_barcode_length:
            argparse.ArgumentError("An i7 index fastq file must be given to attach a sample barcode")
        if args.cell_barcode_length and args.molecule_barcode_length:
            cls._validate_barcode_input(args.molecule_barcode_start_pos,
                                        args.cell_barcode_start_pos + args.cell_barcode_length)

        if args.cell_barcode_length:
            cls.cell_barcode = cls._get_barcode(args.cell_barcode_start_pos,
                                                args.cell_barcode_length,
                                                consts.QUALITY_CELL_BARCODE_TAG_KEY,
                                                consts.RAW_CELL_BARCODE_TAG_KEY)
        if args.molecule_barcode_length:
            cls.molecule_barcode = cls._get_barcode(args.molecule_barcode_start_pos,
                                                args.molecule_barcode_length,
                                                consts.QUALITY_MOLECULE_BARCODE_TAG_KEY,
                                                consts.RAW_MOLECULE_BARCODE_TAG_KEY)
        if args.sample_barcode_length:
            cls.sample_barcode = cls._get_barcode(args.sample_barcode_start_pos,
                                                    args.sample_barcode_length,
                                                    consts.QUALITY_SAMPLE_BARCODE_TAG_KEY,
                                                    consts.RAW_SAMPLE_BARCODE_TAG_KEY)

        tag_generators = cls._make_tag_generators(args.r1, args.i1, args.whitelist)
        cls._tag_bamfile(args.u2, args.output_bamfile, tag_generators)

        return 0
