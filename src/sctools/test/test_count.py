"""
Testing for Count Matrix Construction
=====================================

The test generates (1) a random count matrix, and (2) corresponding alignment records, and writes them to disk
(a BAM file, count matrix, row and column indices). The alignment records are expected to produce the same count
matrix according to the counting algorithm implemented in `sctools:bam.from_sorted_tagged_bam`. Gene names are
fetched from an annotations GTF file that is a subset of GENCODE annotations (see `_test_annotation_file` below).

Notes
-----

- The agreement between the synthetic count matrix and the synthetic BAM file is contingent on the
  agreement between the counting algorithm implemented in `sctools:bam.from_sorted_tagged_bam` and
  the test data generator (see SyntheticTaggedBAMGenerator below). Therefore, future changes in the
  counting algorithm must be accompanied by a corresponding change in the test data generation class.
  Otherwise, the tests will fail.

- We have adopt a minimal test suite design strategy, in the sense that the synthetic test data is only complete
  to the degree that is required by `sctools:bam.from_sorted_tagged_bam`. As such, the synthetic BAM file lacks
  the following features:

    * flag,
    * query_sequence,
    * query_quality,
    * CIGAR string,
    * cell barcode quality tag,
    * molecule barcode quality tag,
    * raw cell and molecule barcodes,

  At the time of writing, the counting algorithm **only** relies on the BAM tags.

- SyntheticTaggedBAMGenerator generates four types of alignment records:

    * necessary alignments -- these records contain one unique cell/molecule/gene tag for each cell/gene count
      unit, according to the randomly generated count matrix. Necessary alignments are also sufficient
      in the sense that they are expected to reproduce the count matrix in the absence of any other alignment
      record.

    * redundant alignments -- these records are expected to be ignored by the counting algorithm and have three
      subtypes:

        - duplicate alignments -- these are randomly picked from necessary alignments, though, they are given a new
          query name (to mimic PCR and optical duplicates).

        - incomplete alignments -- these records miss at least one necessary tag, e.g. cell barcode, molecule
          barcode, or gene name.

        - multi-gene alignments -- these records have the same tags and query_name, though, at least two such
          records per query_name exist that point to different genes.
"""

import operator
import os
import tempfile
from typing import Callable, Optional, List, Set, Tuple, Dict, Generator

import numpy as np
import scipy.sparse as sp
import pysam
import pytest

from sctools import gtf, bam, consts
from sctools.count import CountMatrix

# set the input and output directories
_test_data_dir = os.path.join(os.path.split(__file__)[0], "data")
_test_annotation_file = os.path.join(_test_data_dir, "chr1.30k_records.gtf.gz")

# constants
_test_num_cells = 50
_test_max_genes = 20
_test_gene_expression_rate = 5.0
_test_num_duplicates = 20
_test_num_missing_some_tags = 20
_test_num_multiple_gene_alignments = 20
_test_max_gene_hits_per_multiple_gene_alignments = 5

_test_num_only_exons = 10
_test_num_only_introns = 10
_test_both_exons_introns = 10




@pytest.fixture(scope="module")
def gene_name_to_index() -> Dict[str, int]:
    return gtf.extract_gene_names(_test_annotation_file)


class AlignmentRecordTags:
    """Represents the bundle of cell barcode, molecule barcode, and gene name."""

    def __init__(
        self,
        cell_barcode: Optional[str],
        molecule_barcode: Optional[str],
        gene_name: Optional[str],
        alignment_location: Optional[str] = "EXONIC"
    ) -> None:
        self.cell_barcode = cell_barcode
        self.molecule_barcode = molecule_barcode
        self.gene_name = gene_name
        self.alignment_location = alignment_location

    def __hash__(self):
        return hash((self.cell_barcode, self.molecule_barcode, self.gene_name))

    def __repr__(self):
        return (
            f"{consts.CELL_BARCODE_TAG_KEY}: {self.cell_barcode}, "
            f"{consts.MOLECULE_BARCODE_TAG_KEY}: {self.molecule_barcode}, "
            f"{consts.GENE_NAME_TAG_KEY}: {self.gene_name}",
            f"{consts.ALIGNMENT_LOCATION_TAG_KEY}: {self.alignment_location}"
        )


class CellMoleculeGeneQueryNameSortOrder(bam.AlignmentSortOrder):
    """Hierarchical alignment record sort order (cell barcode >= molecule barcode >= gene name >= query name)."""

    def __init__(
        self,
        cell_barcode_tag_key: str = consts.CELL_BARCODE_TAG_KEY,
        molecule_barcode_tag_key: str = consts.MOLECULE_BARCODE_TAG_KEY,
        gene_name_tag_key: str = consts.GENE_NAME_TAG_KEY,
    ) -> None:
        assert cell_barcode_tag_key, "Cell barcode tag key can not be None"
        assert molecule_barcode_tag_key, "Molecule barcode tag key can not be None"
        assert gene_name_tag_key, "Gene name tag key can not be None"
        self.cell_barcode_tag_key = cell_barcode_tag_key
        self.molecule_barcode_tag_key = molecule_barcode_tag_key
        self.gene_name_tag_key = gene_name_tag_key

    def _get_sort_key(
        self, alignment: pysam.AlignedSegment
    ) -> Tuple[str, str, str, str]:
        return (
            bam.get_tag_or_default(alignment, self.cell_barcode_tag_key, default="N"),
            bam.get_tag_or_default(
                alignment, self.molecule_barcode_tag_key, default="N"
            ),
            bam.get_tag_or_default(alignment, self.gene_name_tag_key, default="N"),
            alignment.query_name,
        )

    @property
    def key_generator(
        self,
    ) -> Callable[[pysam.AlignedSegment], Tuple[str, str, str, str]]:
        return self._get_sort_key

    def __repr__(self) -> str:
        return "hierarchical__cell_molecule_gene_query_name"


class SyntheticTaggedBAMGenerator:
    """This class generates a synthetic count matrix and an accompanying synthetic tagged BAM file as
    described in the preamble documentation block.

    Parameters
    ----------
    num_cells : int
        number of real cells
    max-genes : int
        maximum number of genes to use to generate synthetic counts
    gene_name_to_index : dict
        a map from gene name to their count matrix index
    gene_expression_rate : float
        poisson rate at which each gene is expressed
    rng_seed : int
        random number generator seed

    Methods
    -------
    generate_synthetic_bam_and_counts_matrix
        generates synthetic test data and writes the output to disk

    See Also
    --------
    count.from_sorted_tagged_bam
    """

    OUTPUT_PREFIX = "synthetic_"
    SYNTHETIC_SEQUENCE_NAME = "SYNTHETIC_SEQUENCE"
    SYNTHETIC_SEQUENCE_LENGTH = 100
    NECESSARY_QUERY_NAME_PREFIX = "NECESSARY_QUERY_"
    DUPLICATE_QUERY_NAME_PREFIX = "DUPLICATE_QUERY_"
    INCOMPLETE_QUERY_NAME_PREFIX = "INCOMPLETE_QUERY_"
    MULTI_GENE_QUERY_NAME_PREFIX = "MULTI_GENE_QUERY_"

    bam_output_filename = OUTPUT_PREFIX + "records.bam"
    count_matrix_output_filename = OUTPUT_PREFIX + "count_matrix.npy"
    row_index_output_filename = OUTPUT_PREFIX + "_row_index.npy"
    col_index_output_filename = OUTPUT_PREFIX + "_col_index.npy"

    def __init__(
        self,
        num_cells: int,
        max_genes: int,
        gene_name_to_index: Dict[str, int],
        gene_expression_rate: float,
        rng_seed: int = 777,
    ) -> None:
        self.num_cells = num_cells
        self.gene_expression_rate = gene_expression_rate

        # initialize the random number generator
        self.rng: np.random.RandomState = np.random.RandomState(seed=rng_seed)

        # generate gene names
        self.all_gene_names = [
            k for k, v in sorted(gene_name_to_index.items(), key=operator.itemgetter(1))
        ]
        self.num_genes = len(self.all_gene_names)

        self.max_genes = max_genes
        assert (
            max_genes <= self.num_genes
        ), f"Max genes ({self.max_genes}) must be <= to all annotated genes ({self.num_genes})"
        self.to_be_used_gene_indices: List[int] = self.rng.choice(
            np.arange(0, self.num_genes, dtype=np.int),
            size=self.max_genes,
            replace=False,
        ).tolist()
        self.to_be_used_gene_names = [
            self.all_gene_names[j] for j in self.to_be_used_gene_indices
        ]

    def generate_synthetic_bam_and_counts_matrix(
        self,
        output_path: str,
        num_duplicates: int,
        num_missing_some_tags: int,
        num_multiple_gene_alignments: int,
        max_gene_hits_per_multiple_gene_alignments: int,
        alignment_sort_order: bam.AlignmentSortOrder = CellMoleculeGeneQueryNameSortOrder(),
    ):
        """Generates synthetic count matrix and BAM file and writes them to disk.

        Parameters
        ----------
        output_path : str
            output path
        num_duplicates : int
            number of duplicate records
        num_missing_some_tags : int
            number of records that miss at least one crucial tag
        num_multiple_gene_alignments : int
            number of records that have at least two different gene tags
        max_gene_hits_per_multiple_gene_alignments : int
            maximum number of unique gene names to use for multiple-gene records
        alignment_sort_order : bam.AlignmentSortOrder
            sort order of BAM alignment records; if 'None', random sort order is implied

        Returns
        -------
        None
        """
        assert 2 <= max_gene_hits_per_multiple_gene_alignments <= self.max_genes, (
            f"The parameter `max_gene_hits_per_multiple_gene_alignments` must >= 2 and < maximum annotated "
            f"genes ({self.max_genes})"
        )
        assert num_duplicates >= 0, "Number of duplicate queries must be non-negative"
        assert (
            num_missing_some_tags >= 0
        ), "Number of queries with missing tags must be non-negative"
        assert (
            num_multiple_gene_alignments >= 0
        ), "Number of queries with multiple gene alignments must be non-negative"

        # generate synthetic count matrix and corresponding simulated records
        synthetic_data_bundle = self._generate_synthetic_counts_and_alignment_tags(
            num_duplicates,
            num_missing_some_tags,
            num_multiple_gene_alignments,
            max_gene_hits_per_multiple_gene_alignments,
        )
        records = list(
            SyntheticTaggedBAMGenerator._get_bam_records_generator(
                synthetic_data_bundle
            )
        )

        if not alignment_sort_order:  # random
            # shuffle records
            self.rng.shuffle(records)

        else:
            records = sorted(records, key=alignment_sort_order.key_generator)

        # write BAM file
        with pysam.AlignmentFile(
            os.path.join(output_path, self.bam_output_filename),
            mode="wb",
            reference_names=[self.SYNTHETIC_SEQUENCE_NAME],
            reference_lengths=[self.SYNTHETIC_SEQUENCE_LENGTH],
        ) as bo:
            for record in records:
                bo.write(record)

        # write count matrix, row index, and col index
        np.save(
            os.path.join(output_path, self.count_matrix_output_filename),
            synthetic_data_bundle.count_matrix,
        )
        np.save(
            os.path.join(output_path, self.row_index_output_filename),
            synthetic_data_bundle.row_index,
        )
        np.save(
            os.path.join(output_path, self.col_index_output_filename),
            synthetic_data_bundle.col_index,
        )

    def _generate_synthetic_counts_and_alignment_tags(
        self,
        num_duplicates: int,
        num_missing_some_tags: int,
        num_multiple_gene_alignments: int,
        max_gene_hits_per_multiple_gene_alignments: int,
    ) -> "SyntheticDataBundle":

        # generate count matrix
        count_matrix: np.ndarray = self._generate_random_count_matrix()

        # generate necessary alignment tags that produce count_matrix
        (
            necessary_alignment_record_tags_set,
            row_index,
            col_index,
        ) = self._generate_necessary_alignment_record_bundle(count_matrix)
        necessary_alignment_record_tags_list = list(necessary_alignment_record_tags_set)

        # sanity check -- we require as many necessary alignment records as the total counts
        assert len(necessary_alignment_record_tags_set) == np.sum(count_matrix), (
            "There is an inconsistency between synthetic counts and necessary tags: we require as "
            "many necessary alignment tags as the total counts"
        )

        # add duplicate records
        duplicate_alignment_tags_list = self._generate_duplicate_alignment_tags(
            num_duplicates, necessary_alignment_record_tags_list
        )

        # add records with missing tags
        incomplete_alignment_tags_list: List[
            AlignmentRecordTags
        ] = self._generate_incomplete_alignment_tags(num_missing_some_tags)

        # add records with multiple gene alignments
        multiple_alignment_tags_list: List[
            List[AlignmentRecordTags]
        ] = self._generate_multiple_gene_alignment_tags(
            num_multiple_gene_alignments,
            max_gene_hits_per_multiple_gene_alignments,
            necessary_alignment_record_tags_set,
        )

        return SyntheticDataBundle(
            count_matrix,
            row_index,
            col_index,
            necessary_alignment_record_tags_list,
            duplicate_alignment_tags_list,
            incomplete_alignment_tags_list,
            multiple_alignment_tags_list,
        )

    def _generate_random_count_matrix(self) -> np.ndarray:
        """Generates a random count matrix.

        This method selects `self.max_genes` out of all all genes (`self.num_genes`) and populates the selected genes
        with Poisson counts with rate `self.gene_expression_rate`. The count matrix entries corresponding to the
        rest of the genes are set to zero.

        Returns
        -------
        np.ndarray
            an ndarray of shape (`self.num_cells`, `self.num_genes`)
        """
        non_zero_count_matrix = self.rng.poisson(
            lam=self.gene_expression_rate, size=(self.num_cells, self.max_genes)
        )
        count_matrix = np.zeros((self.num_cells, self.num_genes), dtype=np.int)
        for i, i_gene in enumerate(self.to_be_used_gene_indices):
            count_matrix[:, i_gene] = non_zero_count_matrix[:, i]
        return count_matrix

    @staticmethod
    def _get_bam_records_generator(
        synthetic_data_bundle: "SyntheticDataBundle", rng_seed: int = 777
    ) -> Generator[pysam.AlignedSegment, None, None]:
        """Returns a generator of pysam.AlignedSegment instances created from the alignment tags
        provided to the initializer.

        Parameters
        ----------
        synthetic_data_bundle : SyntheticDataBundle
            a bundle of synthetic alignment tags
        rng_seed : int
            random number generator seed; it is used for generating random reference_start position.

        See Also
        --------
        - The preamble documentation block for a description of the meaning of different alignment records
          (necessary, duplicate, incomplete, etc.)
        - SyntheticTaggedBAMGenerator._generate_aligned_segment_from_tags
        """
        rng = np.random.RandomState(rng_seed)

        num_queries = synthetic_data_bundle.num_queries
        i_query = 0

        # necessary, duplicate, and incomplete alignments
        for alignment_tags_list, query_name_prefix in zip(
            [
                synthetic_data_bundle.necessary_alignment_record_tags_list,
                synthetic_data_bundle.duplicate_alignment_tags_list,
                synthetic_data_bundle.incomplete_alignment_tags_list,
            ],
            [
                SyntheticTaggedBAMGenerator.NECESSARY_QUERY_NAME_PREFIX,
                SyntheticTaggedBAMGenerator.DUPLICATE_QUERY_NAME_PREFIX,
                SyntheticTaggedBAMGenerator.INCOMPLETE_QUERY_NAME_PREFIX,
            ],
        ):
            for alignment_tags in alignment_tags_list:
                yield SyntheticTaggedBAMGenerator._generate_aligned_segment_from_tags(
                    alignment_tags, query_name_prefix, i_query, num_queries, rng
                )
                i_query += 1

        # multi-gene alignments
        for alignment_tags_list in synthetic_data_bundle.multiple_alignment_tags_list:
            # multiple alignments have the same query name (by definition)
            for alignment_tags in alignment_tags_list:
                yield SyntheticTaggedBAMGenerator._generate_aligned_segment_from_tags(
                    alignment_tags,
                    SyntheticTaggedBAMGenerator.MULTI_GENE_QUERY_NAME_PREFIX,
                    i_query,
                    num_queries,
                    rng,
                )
            i_query += 1

    @staticmethod
    def _generate_aligned_segment_from_tags(
        alignment_tags: AlignmentRecordTags,
        query_prefix: str,
        i_query: int,
        num_queries: int,
        rng: np.random.RandomState,
        record_reference_id: Optional[int] = 0,
        reference_start: Optional[int] = -1
    ) -> pysam.AlignedSegment:
        """Generates pysam.AlignedSegment instances from alignment_tags.

        Parameters
        ----------
        alignment_tags : AlignmentRecordTags
            tags to attach to the instantiated pysam.AlignedSegment
        query_prefix : str
            prefix to use for query name
        i_query : int
            query index
        num_queries: int
            maximum number of queries (only used for pretty-printing the query index)
        rng: np.random.RandomState
            a random number generator

        Notes
        -----
        The query_sequence and query_quality are both empty as these query features are not used for generating
        the counts matrix. Likewise, the flag is currently unset. In the future, once we add a filtering
        policy based on BAM record flags (such as duplicates), this method must be updated accordingly.

        Returns
        -------
        pysam.AlignedSegment
            an instance of pysam.AlignedSegment

        """
        tags = []
        if alignment_tags.cell_barcode:
            tags.append((consts.CELL_BARCODE_TAG_KEY, alignment_tags.cell_barcode, "Z"))
        if alignment_tags.molecule_barcode:
            tags.append(
                (consts.MOLECULE_BARCODE_TAG_KEY, alignment_tags.molecule_barcode, "Z")
            )
        if alignment_tags.gene_name:
            tags.append((consts.GENE_NAME_TAG_KEY, alignment_tags.gene_name, "Z"))

        if alignment_tags.alignment_location:
            tags.append((consts.ALIGNMENT_LOCATION_TAG_KEY, alignment_tags.alignment_location, "Z"))

        record = pysam.AlignedSegment()
        record.query_name = SyntheticTaggedBAMGenerator._generate_query_name(
            query_prefix, i_query, num_queries
        )

        if reference_start == -1:
            record.reference_start = rng.randint(
              low=0, high=SyntheticTaggedBAMGenerator.SYNTHETIC_SEQUENCE_LENGTH
            )
        else:
            record.reference_start = reference_start

        record.reference_id = record_reference_id  # note: we only use one synthetic sequence
        if len(tags) > 0:
            record.set_tags(tags)
        return record

    @staticmethod
    def _generate_query_name(query_prefix: str, i_query: int, num_queries: int) -> str:
        """Returns query name string from query index. We zero-pad the string representation of query
        indices merely for pretty-printing, e.g. 0000, 0001, ..., 9999."""
        num_digits = len(str(num_queries - 1))
        return query_prefix + str(i_query).zfill(num_digits)

    def _generate_necessary_alignment_record_bundle(
        self, count_matrix: np.ndarray
    ) -> Tuple[Set[AlignmentRecordTags], List[str], List[str]]:
        alignments: Set[AlignmentRecordTags] = set()
        used_cell_barcodes: Set[str] = set()

        row_index: List[str] = []
        col_index = self.all_gene_names

        for i_cell in range(self.num_cells):
            # generate a unique cell barcode
            while True:
                cell_barcode = self._generate_random_cell_barcode()
                if cell_barcode not in used_cell_barcodes:
                    break
            row_index.append(cell_barcode)

            for i_gene in self.to_be_used_gene_indices:
                for i_molecule in range(count_matrix[i_cell, i_gene]):
                    # generate a unique alignment tag
                    unique_alignment_tag = self._generate_unique_random_alignment_tag(
                        alignments,
                        gene_name=self.all_gene_names[i_gene],
                        cell_barcode=cell_barcode,
                    )
                    alignments.add(unique_alignment_tag)

        return alignments, row_index, col_index

    def _generate_unique_random_alignment_tag(
        self,
        existing_alignment_tags: Set[AlignmentRecordTags],
        gene_name: str,
        cell_barcode: Optional[str] = None,
        molecule_barcode: Optional[str] = None,
    ) -> AlignmentRecordTags:
        assert (
            gene_name in self.to_be_used_gene_names
        ), f"{gene_name} is not an allowed gene for generating synthetic data"

        while True:
            alignment = AlignmentRecordTags(
                cell_barcode=cell_barcode
                if cell_barcode
                else self._generate_random_cell_barcode(),
                molecule_barcode=molecule_barcode
                if molecule_barcode
                else self._generate_random_molecule_barcode(),
                gene_name=gene_name,
            )
            if alignment not in existing_alignment_tags:
                return alignment

    def _generate_duplicate_alignment_tags(
        self, num_duplicates: int, necessary_alignments_list: List[AlignmentRecordTags]
    ) -> List[AlignmentRecordTags]:
        return self.rng.choice(necessary_alignments_list, size=num_duplicates).tolist()

    def _generate_incomplete_alignment_tags(
        self, num_missing_some_tags: int
    ) -> List[AlignmentRecordTags]:
        """Generates alignments with missing crucial tags.

        Notes
        -----
        This method requires each combination of missing tags to occur at least once and may therefore return lists
        that are longer than `num_missing_some_tags`.
        """
        incomplete_alignment_tags_list: List[AlignmentRecordTags] = list()
        tag_mask_occurrences: Set[int] = set()
        i_entries = 0
        while i_entries < num_missing_some_tags or len(tag_mask_occurrences) < 7:
            tag_mask = self.rng.randint(low=0, high=7)
            tag_mask_occurrences.add(tag_mask)
            gene_name = self.rng.choice(self.to_be_used_gene_names)
            alignment = self._generate_unique_random_alignment_tag(set(), gene_name)
            if not tag_mask & 1:
                alignment.cell_barcode = None
            if not tag_mask & 2:
                alignment.molecule_barcode = None
            if not tag_mask & 4:
                alignment.gene_name = None
            incomplete_alignment_tags_list.append(alignment)
            i_entries += 1
        return incomplete_alignment_tags_list

    def _generate_multiple_gene_alignment_tags(
        self,
        num_multiple_gene_alignments: int,
        max_gene_hits_per_multiple_gene_alignments: int,
        necessary_alignment_record_tags_set: Set[AlignmentRecordTags],
    ) -> List[List[AlignmentRecordTags]]:

        necessary_alignment_record_tags_list = list(necessary_alignment_record_tags_set)

        multiple_gene_alignment_tags_list: List[List[AlignmentRecordTags]] = list()
        for _ in range(num_multiple_gene_alignments):
            random_necessary_alignment = self.rng.choice(
                necessary_alignment_record_tags_list
            )
            random_necessary_cell_barcode: str = random_necessary_alignment.cell_barcode
            novel_molecule_barcode: str = self._generate_unique_random_alignment_tag(
                necessary_alignment_record_tags_set,
                gene_name=random_necessary_alignment.gene_name,
                cell_barcode=random_necessary_cell_barcode,
            ).molecule_barcode
            num_gene_hits = self.rng.randint(
                low=2, high=max_gene_hits_per_multiple_gene_alignments + 1
            )
            gene_name_hits = self.rng.choice(
                self.to_be_used_gene_names, replace=False, size=num_gene_hits
            )
            multiple_gene_alignment_tags_list.append(
                [
                    AlignmentRecordTags(
                        random_necessary_cell_barcode, novel_molecule_barcode, gene_name
                    )
                    for gene_name in gene_name_hits
                ]
            )
        return multiple_gene_alignment_tags_list

    def _generate_random_cell_barcode(self, length: int = 16):
        return self._generate_random_genomic_sequences(length)

    def _generate_random_molecule_barcode(self, length: int = 10):
        return self._generate_random_genomic_sequences(length)

    def _generate_random_genomic_sequences(self, length: int):
        return "".join(self.rng.choice(["A", "C", "T", "G"], size=length))


class SyntheticDataBundle:
    """A container for synthetic count matrix, row and column indices, and alignment tags.

    Parameters
    ----------
    count_matrix : np.ndarray
        the cell x gene synthetic count matrix
    row_index : List[str]
        list of cell barcodes
    col_index : List[str]
        list of gene names
    necessary_alignment_record_tags_list : List[AlignmentRecordTags]
        list of necessary alignment tags; alignment records made using these tags are expected to produce
        `count_matrix` once processed by the counting algorithm.
    duplicate_alignment_tags_list : List[AlignmentRecordTags]
        list of duplicate alignment tags (a subset of `necessary_alignment_record_tags_list`)
    incomplete_alignment_tags_list : List[AlignmentRecordTags]
        list of incomplete alignment tags (miss at least one of the required tags: cell, molecule, gene)
    multiple_alignment_tags_list : List[List[AlignmentRecordTags]]
        list of lists of multiple alignment tags; each list element is a list of alignment tags with the
        same molecular barcodes, though, with multiple gene names.

    See Also
    --------
    SyntheticBarcodedBAMGenerator
    """

    def __init__(
        self,
        count_matrix: np.ndarray,
        row_index: List[str],
        col_index: List[str],
        necessary_alignment_record_tags_list: List[AlignmentRecordTags],
        duplicate_alignment_tags_list: List[AlignmentRecordTags],
        incomplete_alignment_tags_list: List[AlignmentRecordTags],
        multiple_alignment_tags_list: List[List[AlignmentRecordTags]],
    ) -> None:

        assert count_matrix.shape == (
            len(row_index),
            len(col_index),
        ), "The shape of the count matrix is inconsistent with the provided row/column indices"

        self.count_matrix = count_matrix
        self.row_index = row_index
        self.col_index = col_index

        self.necessary_alignment_record_tags_list = necessary_alignment_record_tags_list
        self.duplicate_alignment_tags_list = duplicate_alignment_tags_list
        self.incomplete_alignment_tags_list = incomplete_alignment_tags_list
        self.multiple_alignment_tags_list = multiple_alignment_tags_list

        self.num_queries = (
            len(necessary_alignment_record_tags_list)
            + len(duplicate_alignment_tags_list)
            + len(incomplete_alignment_tags_list)
            + len(multiple_alignment_tags_list)
        )


def _get_sorted_count_matrix(
    count_matrix: np.ndarray, row_index: np.ndarray, col_index: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sorted the rows and columns of `count_matrix` and the associated row/column indices.

    Parameters
    ----------
    count_matrix : np.ndarray
        a cell x gene count matrix
    row_index : np.ndarray
        row index of the count matrix (i.e. cell barcodes)
    col_index : np.ndarray
        column index of the count matrix (i.e. gene names)

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        row/column sorted count matrix, sorted row index, sorted column index
    """
    sorted_row_indices = [
        idx for idx, _ in sorted(enumerate(row_index), key=operator.itemgetter(1))
    ]
    sorted_col_indices = [
        idx for idx, _ in sorted(enumerate(col_index), key=operator.itemgetter(1))
    ]
    return (
        count_matrix[sorted_row_indices, :][:, sorted_col_indices],
        row_index[sorted_row_indices],
        col_index[sorted_col_indices],
    )

@pytest.mark.parametrize(
    "alignment_sort_order",
    [bam.QueryNameSortOrder(), CellMoleculeGeneQueryNameSortOrder()],
    ids=["query_name_sort_order", "cell_molecule_gene_query_name_sort_order"],
)
def test_count_matrix_from_bam(
    alignment_sort_order: bam.AlignmentSortOrder, gene_name_to_index
):
    # instantiate a test data generator
    synthetic_data_generator = SyntheticTaggedBAMGenerator(
        _test_num_cells, _test_max_genes, gene_name_to_index, _test_gene_expression_rate
    )

    _test_temp_dir = tempfile.TemporaryDirectory()
    try:
        # generate test data
        synthetic_data_generator.generate_synthetic_bam_and_counts_matrix(
            _test_temp_dir.name,
            _test_num_duplicates,
            _test_num_missing_some_tags,
            _test_num_multiple_gene_alignments,
            _test_max_gene_hits_per_multiple_gene_alignments,
            alignment_sort_order=alignment_sort_order,
        )

        # test data paths
        test_bam_path = os.path.join(
            _test_temp_dir.name, SyntheticTaggedBAMGenerator.bam_output_filename
        )
        test_count_matrix_path = os.path.join(
            _test_temp_dir.name,
            SyntheticTaggedBAMGenerator.count_matrix_output_filename,
        )
        test_row_index_path = os.path.join(
            _test_temp_dir.name, SyntheticTaggedBAMGenerator.row_index_output_filename
        )
        test_col_index_path = os.path.join(
            _test_temp_dir.name, SyntheticTaggedBAMGenerator.col_index_output_filename
        )

        # create CountMatrix from the synthetic bam
        count_matrix_from_bam: CountMatrix = CountMatrix.from_sorted_tagged_bam(
            test_bam_path, gene_name_to_index
        )

        # load the test counts matrix
        count_matrix_data_expected = np.load(test_count_matrix_path)
        row_index_expected = np.load(test_row_index_path)
        col_index_expected = np.load(test_col_index_path)

    finally:
        _test_temp_dir.cleanup()

    count_matrix_data_from_bam = count_matrix_from_bam.matrix.todense()
    row_index_from_bam = count_matrix_from_bam.row_index
    col_index_from_bam = count_matrix_from_bam.col_index

    # sort expected and from_bam results by their respective row and column indices, since their sort order
    # is not part of the design specs and is considered arbitrary
    (
        sorted_count_matrix_data_from_bam,
        sorted_row_index_from_bam,
        sorted_col_index_from_bam,
    ) = _get_sorted_count_matrix(
        count_matrix_data_from_bam, row_index_from_bam, col_index_from_bam
    )
    (
        sorted_count_matrix_data_expected,
        sorted_row_index_expected,
        sorted_col_index_expected,
    ) = _get_sorted_count_matrix(
        count_matrix_data_expected, row_index_expected, col_index_expected
    )

    # assert equality of sorted count matrices and sorted row/col indices
    assert np.allclose(
        sorted_count_matrix_data_from_bam, sorted_count_matrix_data_expected
    )
    assert all(
        [
            row_name_from_bam == row_name_expected
            for row_name_from_bam, row_name_expected in zip(
                sorted_row_index_from_bam, sorted_row_index_expected
            )
        ]
    )
    assert all(
        [
            col_name_from_bam == col_name_expected
            for col_name_from_bam, col_name_expected in zip(
                sorted_col_index_from_bam, sorted_col_index_expected
            )
        ]
    )

@pytest.mark.parametrize(
    "alignment_sort_order",
    [bam.QueryNameSortOrder(), CellMoleculeGeneQueryNameSortOrder()],
    ids=["query_name_sort_order", "cell_molecule_gene_query_name_sort_order"],
)
def test_count_matrix_with_introns(
    alignment_sort_order: bam.AlignmentSortOrder, gene_name_to_index
):
    _count_matrix_with_introns(alignment_sort_order, gene_name_to_index, 0)
    _count_matrix_with_introns(alignment_sort_order, gene_name_to_index, 1)

def _count_matrix_with_introns(
    alignment_sort_order: bam.AlignmentSortOrder, gene_name_to_index, test_index
):

    chromosomes_gene_locations_extended = gtf.extract_extended_gene_names(_test_annotation_file)
    chromosomes_gene_exons = gtf.extract_gene_exons(_test_annotation_file)
    _test_chromosomes_gene_non_exons = gtf.extract_gene_non_exons(chromosomes_gene_exons,
                                                                  chromosomes_gene_locations_extended
                                                                  )

    _test_chromosomes_gene_exons = {}
    for chromosome in chromosomes_gene_exons:
        _test_chromosomes_gene_exons[chromosome]={}
        for gene_exons in chromosomes_gene_exons[chromosome]:
            _test_chromosomes_gene_exons[chromosome][gene_exons[1]] = gene_exons[0]

    # instantiate a test data generator
    chromosome = list(_test_chromosomes_gene_exons.keys())[0]

    synthetic_data_generator = SyntheticTaggedAlignmentTypeBAMGenerator(
        _test_num_cells, _test_max_genes, _test_chromosomes_gene_exons[chromosome],
        _test_chromosomes_gene_non_exons[chromosome]
    )

    _test_temp_dir = tempfile.TemporaryDirectory()
    try:
        # generate test data
        filename= synthetic_data_generator.generate_synthetic_bam_and_counts_matrix(
            _test_temp_dir.name,
            _test_num_only_exons,
            _test_num_only_introns,
            _test_both_exons_introns,
            gene_name_to_index,
            test_index,
            alignment_sort_order=alignment_sort_order
        )

        # test data paths
        test_bam_path = os.path.join(
            _test_temp_dir.name, SyntheticTaggedAlignmentTypeBAMGenerator.bam_output_filename
        )
        test_count_matrix_path = os.path.join(
            _test_temp_dir.name,
            SyntheticTaggedAlignmentTypeBAMGenerator.count_matrix_output_filename,
        )
        test_row_index_path = os.path.join(
            _test_temp_dir.name, SyntheticTaggedAlignmentTypeBAMGenerator.row_index_output_filename
        )
        test_col_index_path = os.path.join(
            _test_temp_dir.name, SyntheticTaggedAlignmentTypeBAMGenerator.col_index_output_filename
        )
        # create CountMatrix from the synthetic bam
        if test_index==0:
           count_matrix_from_bam: CountMatrix = CountMatrix.from_sorted_tagged_bam(
              test_bam_path, gene_name_to_index
           )
        if test_index==1:
            count_matrix_from_bam: CountMatrix = CountMatrix.from_sorted_tagged_bam(
               test_bam_path, gene_name_to_index,
               chromosomes_gene_locations_extended=chromosomes_gene_locations_extended
            )

        # load the test counts matrix
        _count_matrix_data_expected = sp.csr_matrix(np.load(test_count_matrix_path))
        row_index_expected = np.load(test_row_index_path)
        col_index_expected = np.load(test_col_index_path)

        count_matrix_data_expected = CountMatrix(_count_matrix_data_expected, row_index_expected, col_index_expected)
        count_matrix_data_expected = count_matrix_data_expected.matrix.todense()

    finally:
        _test_temp_dir.cleanup()

    count_matrix_data_from_bam = count_matrix_from_bam.matrix.todense()
    row_index_from_bam = count_matrix_from_bam.row_index
    col_index_from_bam = count_matrix_from_bam.col_index

    # sort expected and from_bam results by their respective row and column indices, since their sort order
    # is not part of the design specs and is considered arbitrary
    (
        sorted_count_matrix_data_from_bam,
        sorted_row_index_from_bam,
        sorted_col_index_from_bam,
    ) = _get_sorted_count_matrix(
        count_matrix_data_from_bam, row_index_from_bam, col_index_from_bam
    )
    (
        sorted_count_matrix_data_expected,
        sorted_row_index_expected,
        sorted_col_index_expected,
    ) = _get_sorted_count_matrix(
        count_matrix_data_expected, row_index_expected, col_index_expected
    )

    assert all(
        [
            row_name_from_bam == row_name_expected
            for row_name_from_bam, row_name_expected in zip(
            sorted_row_index_from_bam, sorted_row_index_expected
        )
        ]
    )
    assert all(
        [
            col_name_from_bam == col_name_expected
            for col_name_from_bam, col_name_expected in zip(
            sorted_col_index_from_bam, sorted_col_index_expected
        )
        ]
    )

    assert np.allclose(
        sorted_count_matrix_data_from_bam, sorted_count_matrix_data_expected
    )





class SyntheticTaggedAlignmentTypeBAMGenerator:
    """This class generates a synthetic count matrix and an accompanying synthetic tagged BAM file as
    described in the preamble documentation block.

    Parameters
    ----------
    num_cells : int
        number of real cells
    max-genes : int
        maximum number of genes to use to generate synthetic counts



    rng_seed : int
        random number generator seed

    Methods
    -------
    generate_synthetic_bam_and_counts_matrix
        generates synthetic test data and writes the output to disk

    See Also
    --------
    count.from_sorted_tagged_bam
    """

    OUTPUT_PREFIX = "intronic_"
    SYNTHETIC_SEQUENCE_LENGTH = 5
    REFERENCE_SEQUENCE_NAME = "1"
    EXONIC_SEQUENCE_NAME = "EXONIC_SEQUENCE"
    SYNTHETIC_SEQUENCE_LENGTH = 100
    EXONIC_QUERY_NAME_PREFIX = "EXONIC_QUERY_"
    INTRONIC_QUERY_NAME_PREFIX = "INTRONIC_QUERY_"
    EXONIC_INTRONIC_QUERY_NAME_PREFIX = "EXONIC_INTRONIC_QUERY_"

    bam_output_filename = OUTPUT_PREFIX + "records.bam"
    count_matrix_output_filename = OUTPUT_PREFIX + "count_matrix.npy"
    row_index_output_filename = OUTPUT_PREFIX + "_row_index.npy"
    col_index_output_filename = OUTPUT_PREFIX + "_col_index.npy"

    def __init__(
        self,
        num_cells: int,
        max_genes: int,
        chromosomes_gene_exons: Dict[str, Dict[str, List[tuple]]],
        chromosomes_gene_non_exons: Dict[str, List[tuple]],
        rng_seed: int = 777,
    ) -> None:
        self.num_cells = num_cells

        self.chromosomes_gene_exons = chromosomes_gene_exons
        self.chromosomes_gene_non_exons = chromosomes_gene_non_exons

        # initialize the random number generator
        self.rng: np.random.RandomState = np.random.RandomState(seed=rng_seed)

        # generate gene names
        self.all_gene_names = list(self.chromosomes_gene_exons.keys())[:max_genes]
        self.num_genes = len(self.all_gene_names)

        self.max_genes = max_genes
        assert (
            max_genes <= self.num_genes
        ), f"Max genes ({self.max_genes}) must be <= to all annotated genes ({self.num_genes})"
        self.to_be_used_gene_indices: List[int] = self.rng.choice(
            np.arange(0, self.num_genes, dtype=np.int),
            size=self.max_genes,
            replace=False,
        ).tolist()
        self.to_be_used_gene_names = [
            self.all_gene_names[j] for j in self.to_be_used_gene_indices
        ]

    def _generate_random_cell_barcode(self, length: int = 16):
        return self._generate_random_genomic_sequences(length)

    def _generate_random_molecule_barcode(self, length: int = 10):
        return self._generate_random_genomic_sequences(length)

    def _generate_random_genomic_sequences(self, length: int):
        return "".join(self.rng.choice(["A", "C", "T", "G"], size=length))

    def _generate_location_based_tag_list(
            self,
            num_alignments: int,
            gene_names: List[str],
            alignment_location: str
    ):
        alignment_record_tags = []
        for i in range(num_alignments):
            alignment_record_tags.append(AlignmentRecordTags(
                                                        self._generate_random_cell_barcode(),
                                                        self._generate_random_molecule_barcode(),
                                                        gene_names[i],
                                                        alignment_location
                                                        )
            )

        return alignment_record_tags

    def _add_alignment_start_coordinates(self,
        alignment_tags, alignment_location
    ):
        _alignment_tags = []

        for alignment_tag in alignment_tags:
           if alignment_location=='EXONIC':
               if alignment_tag.gene_name in self.chromosomes_gene_exons:
                   coord= self.chromosomes_gene_exons[alignment_tag.gene_name]
                   setattr(alignment_tag, 'coordinate', coord[0][0] + 1)
                   _alignment_tags.append(alignment_tag)

           if alignment_location=='INTRONIC':
               if alignment_tag.gene_name in self.chromosomes_gene_non_exons:
                   coord= self.chromosomes_gene_non_exons[alignment_tag.gene_name]
                   if coord:
                     setattr(alignment_tag, 'coordinate', coord[0][0] +1)
                     alignment_tag.gene_name = ""
                     _alignment_tags.append(alignment_tag)

        return _alignment_tags



    def generate_synthetic_bam_and_counts_matrix(
                self,
                output_path: str,
                num_only_exons: int,
                num_only_introns: int,
                num_both_exons_introns: int,
                gene_name_to_index: int,
                test_index: int,
                alignment_sort_order: bam.AlignmentSortOrder = CellMoleculeGeneQueryNameSortOrder()
        ):
            """Generates synthetic count matrix and BAM file and writes them to disk.

            Parameters
            ----------
            output_path : str
                output path

            alignment_sort_order : bam.AlignmentSortOrder
                sort order of BAM alignment records; if 'None', random sort order is implied

            Returns
            -------
            None
            """
            assert(num_only_exons >= 0), "Number of exons only alignments must be non-negative"
            assert(num_only_introns >= 0), "Number of introns only must be non-negative"
            assert(num_both_exons_introns >= 0), "Number of queries alignment in both introns and exons must be non-negative"

            gene_names_alignments =[]

            for gene_name in sorted(self.chromosomes_gene_non_exons.keys()):
                if self.chromosomes_gene_non_exons[gene_name]:
                    gene_names_alignments.append(gene_name)

            gene_names: List[int] =[]
            cell_ids: List[int] = []

            records = []
            "Only exons, expected in both single-cell and single-nuclei modes"
            exonic_alignment_tags = self._generate_location_based_tag_list(10, gene_names_alignments[0:],'EXONIC')
            exonic_alignment_tags = self._add_alignment_start_coordinates(exonic_alignment_tags,'EXONIC')

            for i, alignment_tag in enumerate(exonic_alignment_tags):
                pysam_alignment = SyntheticTaggedBAMGenerator._generate_aligned_segment_from_tags(
                        alignment_tag,
                        'EXONIC',
                        i,
                        10,
                        self.rng,
                        reference_start=alignment_tag.coordinate
                )
                records.append(pysam_alignment)
                gene_names.append(alignment_tag.gene_name)
                cell_ids.append(alignment_tag.cell_barcode)

            "Only introns only in single-nuclei mode"
            intronic_alignment_tags = self._generate_location_based_tag_list(3, gene_names_alignments[10:],'INTRONIC')
            intronic_alignment_tags = self._add_alignment_start_coordinates(intronic_alignment_tags,'INTRONIC')
            for i, alignment_tag in enumerate(intronic_alignment_tags):
                pysam_alignment = SyntheticTaggedBAMGenerator._generate_aligned_segment_from_tags(
                        alignment_tag,
                        'INTRONIC',
                        i+10,
                        10,
                        self.rng,
                        reference_start=alignment_tag.coordinate
                )
                records.append(pysam_alignment)
                if test_index==1:
                    gene_names.append(gene_names_alignments[i+10])
                    cell_ids.append(alignment_tag.cell_barcode)

            "both intron and exons from the same gene in bost single-cell and single-nuclei modes"
            exonic_alignment_tags = self._generate_location_based_tag_list(10, gene_names_alignments[20:],'EXONIC')
            exonic_alignment_tags = self._add_alignment_start_coordinates(exonic_alignment_tags,'EXONIC')

            _intronic_alignment_tags = self._generate_location_based_tag_list(10, gene_names_alignments[20:],'INTRONIC')
            intronic_alignment_tags = []
            for intronic_tag, exonic_tag in zip(_intronic_alignment_tags, exonic_alignment_tags):
                intronic_tag.cell_barcode= exonic_tag.cell_barcode
                intronic_alignment_tags.append(intronic_tag)
            intronic_alignment_tags = self._add_alignment_start_coordinates(intronic_alignment_tags,'INTRONIC')

            for i, (exonic_alignment_tag, intronic_alignment_tag) in enumerate(zip(exonic_alignment_tags, intronic_alignment_tags)):
                pysam_alignment = SyntheticTaggedBAMGenerator._generate_aligned_segment_from_tags(
                        exonic_alignment_tag,
                        'EXONINTRONSAME',
                        i+20,
                        10,
                        self.rng,
                        reference_start=exonic_alignment_tag.coordinate
                )
                records.append(pysam_alignment)

                pysam_alignment = SyntheticTaggedBAMGenerator._generate_aligned_segment_from_tags(
                        intronic_alignment_tag,
                        'EXONINTRONSAME',
                        i + 20,
                        10,
                        self.rng,
                        reference_start=intronic_alignment_tag.coordinate
                )
                records.append(pysam_alignment)
                cell_ids.append(exonic_alignment_tag.cell_barcode)
                gene_names.append(exonic_alignment_tag.gene_name)


            "both intron and exons from separate genes should not appear in single-cell mode"
            exonic_alignment_tags = self._generate_location_based_tag_list(10, gene_names_alignments[30:],'EXONIC')
            exonic_alignment_tags = self._add_alignment_start_coordinates(exonic_alignment_tags, 'EXONIC')

            _intronic_alignment_tags = self._generate_location_based_tag_list(10, gene_names_alignments[31:],'INTRONIC')
            intronic_alignment_tags = []
            for intronic_tag, exonic_tag in zip(_intronic_alignment_tags, exonic_alignment_tags):
                intronic_tag.cell_barcode= exonic_tag.cell_barcode
                intronic_alignment_tags.append(intronic_tag)
            intronic_alignment_tags = self._add_alignment_start_coordinates(intronic_alignment_tags, 'INTRONIC')

            for i, (exonic_alignment_tag, intronic_alignment_tag) in enumerate(zip(exonic_alignment_tags,  intronic_alignment_tags)):
                pysam_alignment = SyntheticTaggedBAMGenerator._generate_aligned_segment_from_tags(
                    exonic_alignment_tag, 'EXONINTRONSEP',
                        i + 30,
                        10,
                        self.rng,
                        reference_start=exonic_alignment_tag.coordinate
                )
                records.append(pysam_alignment)

                pysam_alignment = SyntheticTaggedBAMGenerator._generate_aligned_segment_from_tags(
                    intronic_alignment_tag, 'EXONINTRONSEP',
                    i + 30,
                    10,
                    self.rng,
                    reference_start=intronic_alignment_tag.coordinate
                )
                records.append(pysam_alignment)

                if test_index==0:
                    cell_ids.append(exonic_alignment_tag.cell_barcode)
                    gene_names.append(exonic_alignment_tag.gene_name)



            # write BAM file
            with pysam.AlignmentFile(
                    os.path.join(output_path, self.bam_output_filename),
                    mode="wb",
                    reference_names=[self.REFERENCE_SEQUENCE_NAME],
                    reference_lengths=[self.SYNTHETIC_SEQUENCE_LENGTH],
            ) as bo:
                for record in records:
                    bo.write(record)


            n_genes = len(gene_name_to_index)
            n_data= len(cell_ids)
            # write count matrix, row index, and col index
            count_matrix = np.zeros((n_data, n_genes), dtype=np.int32)
            for i, (cell_id, gene_name) in enumerate(zip(cell_ids, gene_names)):
                count_matrix[i][gene_name_to_index[gene_name]] = 1

            test_count_matrix_path = os.path.join(
                output_path,
                SyntheticTaggedAlignmentTypeBAMGenerator.count_matrix_output_filename,
            )
            test_row_index_path = os.path.join(
                output_path, SyntheticTaggedAlignmentTypeBAMGenerator.row_index_output_filename
            )
            test_col_index_path = os.path.join(
                output_path, SyntheticTaggedAlignmentTypeBAMGenerator.col_index_output_filename
            )

            np.save(test_count_matrix_path, count_matrix)
            np.save(test_row_index_path, cell_ids)
            gene_rank = [ (gene, rank) for gene, rank in gene_name_to_index.items()]
            gene_rank.sort( key = lambda x: x[1])
            gene_names = [ x[0] for x in gene_rank ]
            np.save(test_col_index_path, gene_names)

            return os.path.join(output_path, self.bam_output_filename)