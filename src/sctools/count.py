"""
Construct Count Matrices
========================

This module defines methods that enable (optionally) distributed construction of count matrices.
This module outputs coordinate sparse matrices that are converted to CSR matrices prior to delivery
for compact storage, and helper functions to convert this format into other commonly used formats.

Methods
-------
from_sorted_tagged_bam(
    bam_file: str, annotation_file: str, cell_barcode_tag: str = consts.CELL_BARCODE_TAG_KEY,
    molecule_barcode_tag: str=consts.MOLECULE_BARCODE_TAG_KEY,
    gene_name_tag: str=consts.GENE_NAME_TAG_KEY, open_mode: str='rb')
from_mtx(matrix_mtx: str, row_index_file: str, col_index_file: str)

Notes
-----
Memory usage of this module can be roughly approximated by the chunk_size parameter in Optimus.
The memory usage is equal to approximately 6*8 bytes per molecules in the file.
"""

import itertools
import operator
from typing import List, Dict, Tuple, Set, Optional, Generator

import numpy as np
import pysam
import scipy.sparse as sp
from scipy.io import mmread

from sctools import consts, bam


class CountMatrix:
    def __init__(
        self, matrix: sp.csr_matrix, row_index: np.ndarray, col_index: np.ndarray
    ):
        self._matrix = matrix
        self._row_index = row_index
        self._col_index = col_index

    @property
    def matrix(self):
        return self._matrix

    @property
    def row_index(self):
        return self._row_index

    @property
    def col_index(self):
        return self._col_index

    @staticmethod
    def _get_alignments_grouped_by_query_name_generator(
        bam_file: str,
        cell_barcode_tag: str,
        molecule_barcode_tag: str,
        open_mode: str = "rb",
    ) -> Generator[
        Tuple[str, Optional[str], Optional[str], List[pysam.AlignedSegment]], None, None
    ]:
        """Iterates through a query_name-sorted BAM file, groups all alignments with the same query name

        Parameters
        ----------
        bam_file : str
            input bam file marked by cell barcode, molecule barcode, and gene ID tags sorted in that
            order
        cell_barcode_tag : str
            Tag that specifies the cell barcode for each read.
        molecule_barcode_tag : str
            Tag that specifies the molecule barcode for each read.

        Returns
        -------
            a generator for tuples (query_name, cell_barcode, molecule_barcode, alignments)
        """
        with pysam.AlignmentFile(bam_file, mode=open_mode) as bam_records:
            for (query_name, grouper) in itertools.groupby(
                bam_records, key=lambda record: record.query_name
            ):
                alignments: List[pysam.AlignedSegment] = list(grouper)
                cell_barcode: Optional[str] = bam.get_tag_or_default(
                    alignments[0], cell_barcode_tag
                )
                molecule_barcode: Optional[str] = bam.get_tag_or_default(
                    alignments[0], molecule_barcode_tag
                )
                yield query_name, cell_barcode, molecule_barcode, alignments

    """Looks through a list of gene location to fine the one that the given read_start ovelaps

    Parameters
    ----------
    gene_locations: Array
        array with gene start end locations and names
    search_start:
        index of gene to start searching form
    search_end:
        index of gene up to which to search to
    read_start:
        position at which the read starts at

    Returns
    -------
        name of gene with overlap or None if no overlap is found

    """

    @classmethod
    def binary_overlap(cls, gene_locations, search_start, search_end, read_start):
        while search_start <= search_end:
            current_gene_index = int((search_start + search_end) / 2)
            if (
                gene_locations[current_gene_index][0][0]
                < read_start
                < gene_locations[current_gene_index][0][1]
            ):
                return gene_locations[current_gene_index][1]
            elif gene_locations[current_gene_index][0][0] < read_start:
                search_start = current_gene_index + 1
            else:
                search_end = current_gene_index - 1
        return None

    # todo add support for generating a matrix of invalid barcodes
    # todo add support for splitting spliced and unspliced reads
    # todo add support for generating a map of cell barcodes
    # todo add the option for stringent checks on the input (e.g. BAM sort order)
    # todo once the stringent checks are in place, safely move on to the hashset-free implementation
    @classmethod
    def from_sorted_tagged_bam(
        cls,
        bam_file: str,
        gene_name_to_index: Dict[str, int],
        gene_locations: List[tuple] = None,
        cell_barcode_tag: str = consts.CELL_BARCODE_TAG_KEY,
        molecule_barcode_tag: str = consts.MOLECULE_BARCODE_TAG_KEY,
        gene_name_tag: str = consts.GENE_NAME_TAG_KEY,
        open_mode: str = "rb",
    ) -> "CountMatrix":
        """Generate a count matrix from a sorted, tagged bam file

        Notes
        -----
        - Input bam file must be sorted by query name.

        - The sort order of the input BAM file is not strictly checked. If the input BAM file not sorted
          by query_name, the output counts will be wrong without any warnings being issued.

        This method returns counts that correspond to both spliced and unspliced reads.

        Description of the algorithm
        ----------------------------
        The implemented counting strategy is intended to closely match that of CellRanger 2.1.1
        (see the references). The following pseudo-code describes the counting algorithm:

        for each query_name (i.e. unique sequenced read):
            - if only a single alignment exists, _consider_ the read
            - if multiple alignments exist,
                - if a unique gene name is associated to all alignments that have a gene name tag,
                  _consider_ the read; otherwise, the read is useless and neglect it
            - if the read is to be _considered_,
                - if the triple (cell barcode, molecule barcode, gene name) is not encountered before,
                  count it as evidence for a unique transcript; otherwise, consider the read as duplicate
                  and neglect it

        Parameters
        ----------
        bam_file : str
            input bam file marked by cell barcode, molecule barcode, and gene ID tags sorted in that
            order
        gene_locations : List[tuple]
            Location of genes
            (default = None)
        cell_barcode_tag : str, optional
            Tag that specifies the cell barcode for each read. Reads without this tag will be ignored
            (default = consts.CELL_BARCODE_TAG_KEY)
        molecule_barcode_tag : str, optional
            Tag that specifies the molecule barcode for each read. Reads without this tag will be
            ignored (default = consts.MOLECULE_BARCODE_TAG_KEY)
        gene_name_tag
            Tag that specifies the gene name for each read. Reads without this tag will be ignored
            (default = consts.GENE_NAME_TAG_KEY)
        gene_name_to_index : dict
            A map from gene names to their counts matrix column index
        open_mode : {'r', 'rb'}, optional
            indicates that the passed file is a bam file ('rb') or sam file ('r') (default = 'rb').

        Returns
        -------
        count_matrix : CountMatrix
            cells x genes sparse count matrix in compressed sparse row format (cells are compressed)

        Notes
        -----
        All matrices produced by this function called on different BAM chunks that share the same annotation
        file can be concatenated using the scipy sparse vstack function, since by definition, the cell barcodes
        contained in different BAM chunks are mutually exclusive. for example:

        >>> import scipy.sparse as sp
        >>> A = sp.coo_matrix([[1, 2], [3, 4]]).tocsr()
        >>> B = sp.coo_matrix([[5, 6]]).tocsr()
        >>> sp.vstack([A, B]).toarray()
        array([[1, 2],
               [3, 4],
               [5, 6]])

        See Also
        --------
        samtools sort (-t parameter):
            C library that can sort files as required.
            http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS

        TagSortBam.CellSortBam:
            WDL task that accomplishes the sorting necessary for this module.
            https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/TagSortBam.wdl

        Relevant parmalinks to the counting algorithm in CellRanger:
        [1] https://github.com/10XGenomics/cellranger/blob/aba5d379169ff0d4bee60e3d100df35752b90383/mro/stages/counter/
                attach_bcs_and_umis/__init__.py
        [2] https://github.com/10XGenomics/cellranger/blob/aba5d379169ff0d4bee60e3d100df35752b90383/lib/rust/
                annotate_reads/src/main.rs
        """
        # map the gene from reach record to an index in the sparse matrix
        n_genes = len(gene_name_to_index)

        # track which tuples (cell_barcode, molecule_barcode, gene_name) we've encountered so far
        observed_cell_molecule_gene_set: Set[Tuple[str, str, str]] = set()

        # COO sparse matrix entries
        data: List[int] = []
        cell_indices: List[int] = []
        gene_indices: List[int] = []

        # track which cells we've seen, and what the current cell number is
        n_cells = 0
        cell_barcode_to_index: Dict[str, int] = {}

        grouped_records_generator = cls._get_alignments_grouped_by_query_name_generator(
            bam_file, cell_barcode_tag, molecule_barcode_tag, open_mode=open_mode
        )

        for (
            query_name,
            cell_barcode,
            molecule_barcode,
            input_alignments,
        ) in grouped_records_generator:

            # modify alignments to include the gene name to the alignments to INTRONIC regions
            if gene_locations:
                alignments = []
                for alignment in input_alignments:
                    if alignment.has_tag("XF"):
                        aln_type = alignment.get_tag("XF")
                        if (
                            alignment.reference_name
                            and aln_type == "INTRONIC"
                            and alignment.reference_name in gene_locations
                        ):
                            gene_name = cls.binary_overlap(
                                gene_locations[alignment.reference_name],
                                0,
                                len(gene_locations[alignment.reference_name]) - 1,
                                alignment.reference_start,
                            )
                            if gene_name:
                                alignment.set_tag("GE", gene_name)
                    alignments.append(alignment)
            else:
                alignments = input_alignments

            # only keep queries w/ well-formed UMIs
            if cell_barcode is None or molecule_barcode is None:
                continue

            if len(alignments) == 1:
                primary_alignment = alignments[0]
                if primary_alignment.has_tag(gene_name_tag):
                    gene_name = primary_alignment.get_tag(gene_name_tag)
                else:
                    continue  # drop query
            else:  # multi-map
                implicated_gene_names: Set[str] = set()
                for alignment in alignments:
                    if alignment.has_tag(gene_name_tag):
                        implicated_gene_names.add(alignment.get_tag(gene_name_tag))
                if len(implicated_gene_names) == 1:  # only one gene
                    gene_name = implicated_gene_names.__iter__().__next__()
                else:
                    continue  # drop query

            if (
                cell_barcode,
                molecule_barcode,
                gene_name,
            ) in observed_cell_molecule_gene_set:
                continue  # optical/PCR duplicate -> drop query
            else:
                observed_cell_molecule_gene_set.add(
                    (cell_barcode, molecule_barcode, gene_name)
                )

            # find the indices that this molecule should correspond to
            gene_index = gene_name_to_index[gene_name]

            # if we've seen this cell before, get its index, else set it
            try:
                cell_index = cell_barcode_to_index[cell_barcode]
            except KeyError:
                cell_index = n_cells
                cell_barcode_to_index[cell_barcode] = n_cells
                n_cells += 1

            # record the molecule data
            data.append(1)  # one count of this molecule
            cell_indices.append(cell_index)
            gene_indices.append(gene_index)

        # convert into coo_matrix
        coordinate_matrix = sp.coo_matrix(
            (data, (cell_indices, gene_indices)),
            shape=(n_cells, n_genes),
            dtype=np.uint32,
        )

        # convert to a csr sparse matrix and return
        col_index = np.asarray(
            [
                k
                for k, v in sorted(
                    gene_name_to_index.items(), key=operator.itemgetter(1)
                )
            ]
        )
        row_index = np.asarray(
            [
                k
                for k, v in sorted(
                    cell_barcode_to_index.items(), key=operator.itemgetter(1)
                )
            ]
        )

        return cls(coordinate_matrix.tocsr(), row_index, col_index)

    def save(self, prefix: str) -> None:
        sp.save_npz(prefix + ".npz", self._matrix, compressed=True)
        np.save(prefix + "_row_index.npy", self._row_index)
        np.save(prefix + "_col_index.npy", self._col_index)

    @classmethod
    def load(cls, prefix: str) -> "CountMatrix":
        matrix = sp.load_npz(prefix + ".npz")
        row_index = np.load(prefix + "_row_index.npy")
        col_index = np.load(prefix + "_col_index.npy")
        return cls(matrix, row_index, col_index)

    @classmethod
    def merge_matrices(cls, input_prefixes: str) -> "CountMatrix":
        col_indices = [np.load(p + "_col_index.npy") for p in input_prefixes]
        row_indices = [np.load(p + "_row_index.npy") for p in input_prefixes]
        matrices = [sp.load_npz(p + ".npz") for p in input_prefixes]

        matrix: sp.csr_matrix = sp.vstack(matrices, format="csr")
        # todo test that col_indices are all same shape
        col_index = col_indices[0]
        row_index = np.concatenate(row_indices)
        return cls(matrix, row_index, col_index)

    @classmethod
    def from_mtx(
        cls, matrix_mtx: str, row_index_file: str, col_index_file: str
    ) -> "CountMatrix":
        """

        Parameters
        ----------
        matrix_mtx : str
            file containing count matrix in matrix market sparse format
        row_index_file : str
            newline delimited row index file
        col_index_file : str
            newline delimited column index file

        Returns
        -------
        CountMatrix
            instance of class
        """
        matrix: sp.csr_matrix = mmread(matrix_mtx).tocsr()
        with open(row_index_file, "r") as fin:
            row_index = np.array(fin.readlines())
        with open(col_index_file, "r") as fin:
            col_index = np.array(fin.readlines())
        return cls(matrix, row_index, col_index)
