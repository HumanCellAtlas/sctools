"""
Construct Count Matrices
========================

This module defines methods that enable (optionally) distributed construction of count matrices.
This module outputs coordinate sparse matrices that are converted to CSR matrices prior to delivery
for compact storage, and helper functions to convert this format into other commonly used formats.

Methods
-------
bam_to_count(bam_file, cell_barcode_tag: str='CB', molecule_barcode_tag='UB', gene_id_tag='GE')

Notes
-----
Memory usage of this module can be roughly approximated by the chunk_size parameter in Optimus.
The memory usage is equal to approximately 6*8 bytes per molecules in the file.
"""

from typing import List, Dict, Tuple
import tempfile
import operator

import numpy as np
import scipy.sparse as sp
import pysam
import gffutils

from sctools import gtf


class CountMatrix:

    def __init__(self, matrix: sp.csr_matrix, row_index: np.ndarray, col_index: np.ndarray):
        self._matrix = matrix
        self._row_index = row_index
        self._col_index = col_index

    @classmethod
    def from_bam(
            cls,
            bam_file: str,
            annotation_file: str,
            cell_barcode_tag: str='CB',
            molecule_barcode_tag: str='UB',
            gene_id_tag: str='GE',
            open_mode: str='rb',
            ):
        """Generate a count matrix from a sorted, tagged bam file

        Input bam file must be sorted by cell, molecule, and gene (where the gene tag varies fastest).
        This module returns reads that correspond to both spliced and unspliced reads.

        Parameters
        ----------
        bam_file : str
            input bam file marked by cell barcode, molecule barcode, and gene ID tags sorted in that
            order
        cell_barcode_tag : str, optional
            Tag that specifies the cell barcode for each read. Reads without this tag will be ignored
            (default = 'CB')
        molecule_barcode_tag : str, optional
            Tag that specifies the molecule barcode for each read. Reads without this tag will be
            ignored (default = 'UB')
        gene_id_tag
            Tag that specifies the gene for each read. Reads without this tag will be ignored
            (default = 'GE')
        annotation_file : str
            gtf annotation file that was used to create gene ID tags. Used to map genes to indices
        open_mode : {'r', 'rb'}, optional
            indicates that the passed file is a bam file ('rb') or sam file ('r') (default = 'rb').

        Returns
        -------
        count_matrix : CountMatrix
            cells x genes sparse count matrix in compressed sparse row format (cells are compressed)

        Notes
        -----
        Any matrices produced by this function that share the same annotation file can be concatenated
        using the scipy sparse vstack function, for example:

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

        """

        # create input arrays
        data: List[int] = []
        cell_indices: List[int] = []
        gene_indices: List[int] = []

        gene_id_to_index: Dict[str, int] = {}
        gtf_reader = gtf.Reader(annotation_file)

        # map the gene from reach record to an index in the sparse matrix
        for gene_index, record in enumerate(gtf_reader.filter(retain_types=['gene'])):
            gene_id = record.get_attribute('gene_id')
            if gene_id is None:
                raise ValueError(
                    'malformed GTF file detected. Record is of type gene but does not have a '
                    '"gene_id" field: %s' % repr(record))
            gene_id_to_index[gene_id] = gene_index

        # track which cells we've seen, and what the current cell number is
        n_cells = 0
        cell_id_to_index: Dict[str, int] = {}

        # process the data
        current_molecule: Tuple[str, str, str] = tuple()

        with pysam.AlignmentFile(bam_file, mode=open_mode) as f:

            for sam_record in f:

                # get the tags that define the record's molecular identity
                try:
                    gene: str = sam_record.get_tag(gene_id_tag)
                    cell: str = sam_record.get_tag(cell_barcode_tag)
                    molecule: str = sam_record.get_tag(molecule_barcode_tag)
                except KeyError:  # if a record is missing any of these, just drop it.
                    continue

                # each molecule is counted only once
                if current_molecule == (gene, cell, molecule):
                    continue

                # find the indices that this molecule should correspond to
                gene_index = gene_id_to_index[gene]

                # if we've seen this cell before, get its index, else set it
                try:
                    cell_index = cell_id_to_index[cell]
                except KeyError:
                    cell_index = n_cells
                    cell_id_to_index[cell] = n_cells
                    n_cells += 1

                # record the molecule data
                data.append(1)  # one count of this molecule
                cell_indices.append(cell_index)
                gene_indices.append(gene_index)

                # set the current molecule
                current_molecule = (gene, cell, molecule)

        # get shape
        gene_number = len(gene_id_to_index)
        cell_number = len(cell_indices)
        shape = (cell_number, gene_number)

        # convert into coo_matrix
        coordinate_matrix = sp.coo_matrix((data, (cell_indices, gene_indices)),
                                          shape=shape, dtype=np.uint32)

        # convert into csr matrix and return
        col_iterable = [k for k, v in sorted(gene_id_to_index.items(), key=operator.itemgetter(1))]
        row_iterable = [k for k, v in sorted(cell_id_to_index.items(), key=operator.itemgetter(1))]
        col_index = np.array(col_iterable)
        row_index = np.array(row_iterable)
        return cls(coordinate_matrix.tocsr(), row_index, col_index)

    # todo add support for generating a matrix of invalid barcodes
    # todo add support for splitting spliced and unspliced reads
    # todo add support for generating a map of cell barcodes

    def save(self, prefix: str):
        sp.save_npz(prefix + '.npz', self._matrix, compressed=True)
        np.save(prefix + '_row_index.npy', self._row_index)
        np.save(prefix + '_col_index.npy', self._col_index)

    @classmethod
    def load(cls, prefix: str):
        matrix = sp.load_npz(prefix + '.npz')
        row_index = np.load(prefix + '_row_index.npy')
        col_index = np.load(prefix + '_col_index.npy')
        return cls(matrix, row_index, col_index)

    def filter_chromosome(self, chromosome: str, gtf_: str):
        dir_: str = tempfile.mkdtemp()
        db = gffutils.create_db(
            gtf_,
            dbfn=dir_ + 'gtf.db',
            merge_strategy='merge',
            force=True,
            keep_order=True,
            sort_attribute_values=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True)

        keep_genes_mask = np.zeros(self._col_index.shape[0])

        # map gene names in chromosome to a set
        genes_on_chromosome = set()
        for gene in db.features_of_type('gene'):
            if gene.chrom == chromosome:
                genes_on_chromosome.add(gene)

        for i, gene in enumerate(self._col_index):
            if gene in genes_on_chromosome:
                keep_genes_mask[i] = 1

        keep_genes_index = np.where(keep_genes_mask[0])
        csc = self._matrix.tocsc()
        arr_subset = csc[:, keep_genes_index]
        col_subset = self._col_index[keep_genes_index]

        return CountMatrix(arr_subset, self._row_index, col_subset)

    @classmethod
    def from_10x_count_matrix(cls):
        raise NotImplementedError
