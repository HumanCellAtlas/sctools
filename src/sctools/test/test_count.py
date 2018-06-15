"""
Testing for Count Matrix Construction
=====================================
"""

import operator
import os
import tempfile
from typing import Optional, List, Set, Tuple, Generator

import numpy as np
import pysam
import pytest
from sctools import gtf
from sctools.count import CountMatrix

# set the input and output directories
_test_data_dir = os.path.join(os.path.split(__file__)[0], 'data')
_test_annotation_file = os.path.join(_test_data_dir, 'chr1.gtf.gz')

# constants
_test_num_cells = 100
_test_max_genes = 50
_test_gene_expression_rate = 5.0
_test_num_duplicates = 100
_test_num_missing_some_tags = 100
_test_num_multiple_gene_alignments = 100
_test_max_gene_hits_per_multiple_gene_alignments = 10


class AlignmentRecordTags:
    """Stores cell barcode, molecule barcode, and gene name."""
    def __init__(self,
                 cell_barcode: Optional[str],
                 molecule_barcode: Optional[str],
                 gene_name: Optional[str]):
        self.cell_barcode = cell_barcode
        self.molecule_barcode = molecule_barcode
        self.gene_name = gene_name

    def __hash__(self):
        return hash((self.cell_barcode, self.molecule_barcode, self.gene_name))

    def __repr__(self):
        return f'CB: {self.cell_barcode}, UB: {self.molecule_barcode}, GE: {self.gene_name}'


class SyntheticDataBundle:
    """A container for synthetic count matrix, row/col indices, and alignment tags. Once these alignment
    tags are converted to full alignment records, they are expected to produce the synthetic count matrix.

    Parameters
    ----------
    count_matrix : np.ndarray
        the cell x gene synthetic count matrix
    row_index : List[str]
        list of cell barcodes
    col_index : List[str]
        list of gene names
    primary_alignment_record_tags_list : List[AlignmentRecordTags]
        list of primary alignment tags
    duplicate_alignment_tags_list : List[AlignmentRecordTags]
        list of duplicate alignment tags
    incomplete_alignment_tags_list : List[AlignmentRecordTags]
        list of incomplete alignment tags
    multiple_alignment_tags_list : List[AlignmentRecordTags]
        list of lists of multiple alignment tags

    Methods
    -------
    get_bam_records_generator:
        generates alignment records to be written to a BAM files

    See Also
    --------
    SyntheticBarcodedBAMGenerator

    """

    synthetic_sequence_name = "SYNTHETIC_SEQUENCE"
    synthetic_sequence_length = 100
    _primary_query_name_prefix = "PRIMARY_QUERY_"
    _duplicate_query_name_prefix = "DUPLICATE_QUERY_"
    _incomplete_query_name_prefix = "INCOMPLETE_QUERY_"
    _multi_gene_query_name_prefix = "MULTI_GENE_QUERY_"

    def __init__(self,
                 count_matrix: np.ndarray,
                 row_index: List[str],
                 col_index: List[str],
                 primary_alignment_record_tags_list: List[AlignmentRecordTags],
                 duplicate_alignment_tags_list: List[AlignmentRecordTags],
                 incomplete_alignment_tags_list: List[AlignmentRecordTags],
                 multiple_alignment_tags_list: List[List[AlignmentRecordTags]]):

        assert count_matrix.shape == (len(row_index), len(col_index))

        self.count_matrix = count_matrix
        self.row_index = row_index
        self.col_index = col_index

        self.primary_alignment_record_tags_list = primary_alignment_record_tags_list
        self.duplicate_alignment_tags_list = duplicate_alignment_tags_list
        self.incomplete_alignment_tags_list = incomplete_alignment_tags_list
        self.multiple_alignment_tags_list = multiple_alignment_tags_list

        self.num_queries = (len(primary_alignment_record_tags_list)
                            + len(duplicate_alignment_tags_list)
                            + len(incomplete_alignment_tags_list)
                            + len(multiple_alignment_tags_list))

    def _get_query_name(self, query_prefix: str, i_query: int) -> str:
        num_digits = len(str(self.num_queries - 1))
        return query_prefix + str(i_query).zfill(num_digits)

    def _generate_aligned_segment_from_tags(self,
                                            alignment_tags: AlignmentRecordTags,
                                            reference_start: int,
                                            query_prefix: str,
                                            i_query: int) -> pysam.AlignedSegment:
        """Generates pysam.AlignedSegment instances from the contained tags.

        Parameters
        ----------
        alignment_tags : AlignmentRecordTags
            tags to attach to the instantiated pysam.AlignedSegment
        reference_start : int
            start position on the synthetic reference
        query_prefix : str
            prefix to use for query name
        i_query : int
            query index

        Notes
        -----
        The query_sequence and query_quality are both empty as these features are not used for generating
        the counts matrix. Likewise, the flag is currently unset. In the future, once we add a filtering
        policy based on BAM record flags (such as duplicates), this method must be updated accordingly.

        Returns
        -------
        pysam.AlignedSegment
            an instance of pysam.AlignedSegment

        """
        tags = []
        if alignment_tags.cell_barcode:
            tags.append((SyntheticBarcodedBAMGenerator.cell_barcode_tag,
                         alignment_tags.cell_barcode, 'Z'))
        if alignment_tags.molecule_barcode:
            tags.append((SyntheticBarcodedBAMGenerator.molecule_barcode_tag,
                         alignment_tags.molecule_barcode, 'Z'))
        if alignment_tags.gene_name:
            tags.append((SyntheticBarcodedBAMGenerator.gene_name_tag,
                         alignment_tags.gene_name, 'Z'))
        record = pysam.AlignedSegment()
        record.query_name = self._get_query_name(query_prefix, i_query)
        record.reference_start = reference_start
        record.reference_id = 0  # note: we only use one synthetic sequence
        if len(tags) > 0:
            record.set_tags(tags)
        return record

    def get_bam_records_generator(self, rng_seed: int = 777) -> Generator[pysam.AlignedSegment, None, None]:
        """Returns a generator of pysam.AlignedSegment instances created from the alignment tags
        provided to the initializer."""
        rng = np.random.RandomState(rng_seed)
        i_query = 0

        # primary alignments
        for alignment_tags in self.primary_alignment_record_tags_list:
            reference_start = rng.randint(low=0, high=self.synthetic_sequence_length)
            yield self._generate_aligned_segment_from_tags(
                alignment_tags, reference_start, self._primary_query_name_prefix, i_query)
            i_query += 1

        # duplicate alignments
        for alignment_tags in self.duplicate_alignment_tags_list:
            reference_start = rng.randint(low=0, high=self.synthetic_sequence_length)
            yield self._generate_aligned_segment_from_tags(
                alignment_tags, reference_start, self._duplicate_query_name_prefix, i_query)
            i_query += 1

        # incomplete alignments
        for alignment_tags in self.incomplete_alignment_tags_list:
            reference_start = rng.randint(low=0, high=self.synthetic_sequence_length)
            yield self._generate_aligned_segment_from_tags(
                alignment_tags, reference_start, self._incomplete_query_name_prefix, i_query)
            i_query += 1

        # multi-gene alignments
        for alignment_tags_list in self.multiple_alignment_tags_list:
            # multiple alignments have the same query name (by definition)
            for alignment_tags in alignment_tags_list:
                reference_start = rng.randint(low=0, high=self.synthetic_sequence_length)
                yield self._generate_aligned_segment_from_tags(
                    alignment_tags, reference_start, self._multi_gene_query_name_prefix, i_query)
            i_query += 1


class SyntheticBarcodedBAMGenerator:
    """This class generates a synthetic count matrix and an accompanying synthetic barcoded BAM file that is
    expected to reproduce the count matrix.

    Notes
    -----
    - The correspondence between the synthetic count matrix and the synthetic BAM file is modeled after the
      algorithm implemented in `count.from_sorted_tagged_bam`. Therefore, as `count.from_sorted_tagged_bam` evolves,
      this test data generator needs to be updated accordingly.

    - At the moment, the synthetic BAM file lacks flag, query_sequence, query_quality, barcode quality tags,
     and raw barcode tags since `count.from_sorted_tagged_bam` does not currently utilize such features.

    - In brief, this class generates four types of alignment records:

      1) primary alignments: these records contain one unique cell/molecule/gene tag for each cell/gene count
         unit, in correspondence with a randomly generated count matrix. As such, the primary alignment records
         alone are expected to reproduce the randomly generated count matrix.

      2) duplicate alignments: randomly picked from primary alignments though given a new query name. These
         record correspond to PCR and optical duplicates. `count.from_sorted_tagged_bam` must ignore such reads.

      3) incomplete alignments: these records miss at least one necessary tag. `count.from_sorted_tagged_bam` must
         ignore such reads.

      4) multi-gene alignments: these records have the same tags and query name, though, at least two such
         records exist that point to different genes. `count.from_sorted_tagged_bam` must ignore such reads.

    Parameters
    ----------
    num_cells : int
        number of real cells
    max-genes : int
        maximum number of genes to use to generate synthetic counts
    annotation_file : str
        annotation GTF file to look up gene names from
    gene_expression_rate : float
        poisson rate at which each gene is expressed
    rng_seed : int
        random number generator seed

    Methods
    -------
    generate_synthetic_counts_and_alignment_tags
        generates synthetic test data as described above and returns an instance of SyntheticDataBundle
    generate_synthetic_bam_and_counts_matrix
        generates synthetic test data and writes them to disk

    See Also
    --------
    count.from_sorted_tagged_bam
    """

    cell_barcode_tag = 'CB'
    molecule_barcode_tag = 'UB'
    gene_name_tag = 'GE'
    output_prefix = 'synthetic_'
    bam_output_filename = output_prefix + 'records.bam'
    count_matrix_output_filename = output_prefix + 'count_matrix.npy'
    row_index_output_filename = output_prefix + '_row_index.npy'
    col_index_output_filename = output_prefix + '_col_index.npy'

    def __init__(self,
                 num_cells: int,
                 max_genes: int,
                 annotation_file: str,
                 gene_expression_rate: float,
                 rng_seed: int = 777):
        self.num_cells = num_cells
        self.gene_expression_rate = gene_expression_rate

        # initialize the random number generator
        self.rng: np.random.RandomState = np.random.RandomState(seed=rng_seed)

        # generate gene names
        self.all_gene_names = self._load_gene_names(annotation_file)
        self.max_genes = max_genes
        self.num_genes = len(self.all_gene_names)
        self.to_be_used_gene_indices: List[int] = self.rng.choice(
            np.arange(0, self.num_genes, dtype=np.int), size=self.max_genes, replace=False).tolist()
        self.to_be_used_gene_names = [self.all_gene_names[j] for j in self.to_be_used_gene_indices]

    def generate_synthetic_counts_and_alignment_tags(self,
                                                     num_duplicates: int,
                                                     num_missing_some_tags: int,
                                                     num_multiple_gene_alignments: int,
                                                     max_gene_hits_per_multiple_gene_alignments: int)\
            -> SyntheticDataBundle:

        # generate count matrix
        count_matrix: np.ndarray = self._generate_random_count_matrix()

        # generate primary alignment tags that produce count_matrix
        primary_alignment_record_tags_set, row_index, col_index = \
            self._generate_primary_alignment_record_bundle(count_matrix)
        primary_alignment_record_tags_list = list(primary_alignment_record_tags_set)

        # sanity check -- we require as many primary alignment records as the total counts
        assert len(primary_alignment_record_tags_set) == np.sum(count_matrix)

        # add duplicate records
        duplicate_alignment_tags_list = self._generate_duplicate_alignment_tags(
            num_duplicates, primary_alignment_record_tags_list)

        # add records with missing tags
        incomplete_alignment_tags_list: List[AlignmentRecordTags] = self._generate_incomplete_alignment_tags(
            num_missing_some_tags)

        # add records with multiple gene alignments
        multiple_alignment_tags_list: List[List[AlignmentRecordTags]] = self._generate_multiple_gene_alignment_tags(
            num_multiple_gene_alignments,
            max_gene_hits_per_multiple_gene_alignments,
            primary_alignment_record_tags_set)

        return SyntheticDataBundle(count_matrix,
                                   row_index,
                                   col_index,
                                   primary_alignment_record_tags_list,
                                   duplicate_alignment_tags_list,
                                   incomplete_alignment_tags_list,
                                   multiple_alignment_tags_list)

    def generate_synthetic_bam_and_counts_matrix(self,
                                                 output_path: str,
                                                 num_duplicates: int,
                                                 num_missing_some_tags: int,
                                                 num_multiple_gene_alignments: int,
                                                 max_gene_hits_per_multiple_gene_alignments: int,
                                                 query_sort_order: str = 'CB_UB_GE_n'):
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
        query_sort_order : str
            'unsorted' : shuffle the records before saving to BAM
            'CB_UB_GE_n' : sort hierarchically by CB, UB, GE, and query_name
            'n' : only sort by query_name

        Returns
        -------
        None
        """
        assert 2 <= max_gene_hits_per_multiple_gene_alignments < self.max_genes
        assert num_duplicates >= 0
        assert num_missing_some_tags >= 0
        assert num_multiple_gene_alignments >= 0

        # generate synthetic count matrix and corresponding simulated records
        synthetic_data_bundle = self.generate_synthetic_counts_and_alignment_tags(
            num_duplicates,
            num_missing_some_tags,
            num_multiple_gene_alignments,
            max_gene_hits_per_multiple_gene_alignments)
        records = [record for record in synthetic_data_bundle.get_bam_records_generator()]

        if query_sort_order == 'random':
            # shuffle records
            self.rng.shuffle(records)

        elif query_sort_order == 'CB_UB_GE_n':
            # sort by (CB, GE, UB, query_name)
            def record_sort_key(record: pysam.AlignedSegment):
                query_name = record.query_name
                try:
                    cell_barcode = record.get_tag(self.cell_barcode_tag)
                except KeyError:
                    cell_barcode = 'N'
                try:
                    molecule_barcode = record.get_tag(self.molecule_barcode_tag)
                except KeyError:
                    molecule_barcode = 'N'
                try:
                    gene_name = record.get_tag(self.gene_name_tag)
                except KeyError:
                    gene_name = 'N'
                return cell_barcode, molecule_barcode, gene_name, query_name
            records = sorted(records, key=record_sort_key)

        elif query_sort_order == 'n':
            def record_sort_key(record: pysam.AlignedSegment):
                return record.query_name
            records = sorted(records, key=record_sort_key)

        else:
            raise Exception(f'Unknown query sort order {query_sort_order}')

        # write BAM file
        with pysam.AlignmentFile(os.path.join(output_path, self.bam_output_filename),
                                 mode='wb',
                                 reference_names=[SyntheticDataBundle.synthetic_sequence_name],
                                 reference_lengths=[SyntheticDataBundle.synthetic_sequence_length]) as bo:
            for record in records:
                bo.write(record)

        # write count matrix, row index, and col index
        np.save(os.path.join(output_path, self.count_matrix_output_filename), synthetic_data_bundle.count_matrix)
        np.save(os.path.join(output_path, self.row_index_output_filename), synthetic_data_bundle.row_index)
        np.save(os.path.join(output_path, self.col_index_output_filename), synthetic_data_bundle.col_index)

    def _generate_random_count_matrix(self) -> np.ndarray:
        non_zero_count_matrix = self.rng.poisson(
            lam=self.gene_expression_rate, size=(self.num_cells, self.max_genes))
        count_matrix = np.zeros((self.num_cells, self.num_genes), dtype=np.int)
        for i, i_gene in enumerate(self.to_be_used_gene_indices):
            count_matrix[:, i_gene] = non_zero_count_matrix[:, i]
        return count_matrix

    def _generate_primary_alignment_record_bundle(self, count_matrix: np.ndarray) \
            -> Tuple[Set[AlignmentRecordTags], List[str], List[str]]:
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
                        cell_barcode=cell_barcode)
                    alignments.add(unique_alignment_tag)

        return alignments, row_index, col_index

    def _generate_unique_random_alignment_tag(self,
                                              existing_alignment_tags: Set[AlignmentRecordTags],
                                              gene_name: str,
                                              cell_barcode: Optional[str] = None,
                                              molecule_barcode: Optional[str] = None) -> AlignmentRecordTags:
        assert gene_name in self.to_be_used_gene_names, f"{gene_name} is not an allowed gene for generating " \
                                                        f"synthetic data"
        while True:
            alignment = AlignmentRecordTags(
                cell_barcode=cell_barcode if cell_barcode else self._generate_random_cell_barcode(),
                molecule_barcode=molecule_barcode if molecule_barcode else self._generate_random_molecule_barcode(),
                gene_name=gene_name)
            if alignment not in existing_alignment_tags:
                return alignment

    def _generate_duplicate_alignment_tags(self,
                                           num_duplicates: int,
                                           primary_alignments_list: List[AlignmentRecordTags])\
            -> List[AlignmentRecordTags]:
        return self.rng.choice(primary_alignments_list, size=num_duplicates).tolist()

    def _generate_incomplete_alignment_tags(self, num_missing_some_tags: int) -> List[AlignmentRecordTags]:
        incomplete_alignment_tags_list: List[AlignmentRecordTags] = list()
        for _ in range(num_missing_some_tags):
            tag_mask = self.rng.randint(low=0, high=7)
            gene_name = self.rng.choice(self.to_be_used_gene_names)
            alignment = self._generate_unique_random_alignment_tag(set(), gene_name)
            if not tag_mask & 1:
                alignment.cell_barcode = None
            if not tag_mask & 2:
                alignment.molecule_barcode = None
            if not tag_mask & 4:
                alignment.gene_name = None
            incomplete_alignment_tags_list.append(alignment)
        return incomplete_alignment_tags_list

    def _generate_multiple_gene_alignment_tags(self,
                                               num_multiple_gene_alignments: int,
                                               max_gene_hits_per_multiple_gene_alignments: int,
                                               primary_alignment_record_tags_set: Set[AlignmentRecordTags]) -> \
            List[List[AlignmentRecordTags]]:

        primary_alignment_record_tags_list = list(primary_alignment_record_tags_set)

        multiple_gene_alignment_tags_list: List[List[AlignmentRecordTags]] = list()
        for _ in range(num_multiple_gene_alignments):
            random_primary_alignment = self.rng.choice(primary_alignment_record_tags_list)
            random_primary_cell_barcode: str = random_primary_alignment.cell_barcode
            novel_molecule_barcode: str = self._generate_unique_random_alignment_tag(
                primary_alignment_record_tags_set,
                gene_name=random_primary_alignment.gene_name,
                cell_barcode=random_primary_cell_barcode).molecule_barcode
            num_gene_hits = self.rng.randint(low=2, high=max_gene_hits_per_multiple_gene_alignments + 1)
            gene_name_hits = self.rng.choice(self.to_be_used_gene_names, replace=False, size=num_gene_hits)
            multiple_gene_alignment_tags_list.append(
                [AlignmentRecordTags(random_primary_cell_barcode, novel_molecule_barcode, gene_name)
                 for gene_name in gene_name_hits])
        return multiple_gene_alignment_tags_list

    @staticmethod
    def _load_gene_names(annotation_file: str) -> List[str]:
        return [k for k, v in sorted(gtf.extract_gene_names(annotation_file).items(), key=operator.itemgetter(1))]

    def _generate_random_cell_barcode(self, length: int = 16):
        return self._generate_random_genomic_sequences(length)

    def _generate_random_molecule_barcode(self, length: int = 10):
        return self._generate_random_genomic_sequences(length)

    def _generate_random_genomic_sequences(self, length: int):
        return ''.join(self.rng.choice(['A', 'C', 'T', 'G'], size=length))


def _get_sorted_count_matrix(count_matrix: np.ndarray, row_index: np.ndarray, col_index: np.ndarray)\
        -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    sorted_row_indices = [idx for idx, _ in sorted(enumerate(row_index), key=operator.itemgetter(1))]
    sorted_col_indices = [idx for idx, _ in sorted(enumerate(col_index), key=operator.itemgetter(1))]
    return (count_matrix[sorted_row_indices, :][:, sorted_col_indices],
            row_index[sorted_row_indices],
            col_index[sorted_col_indices])


@pytest.mark.parametrize('query_sort_order', ['CB_UB_GE_n', 'n'])
def test_count_matrix_from_bam(query_sort_order):
    # instantiate a test data generator
    synthetic_data_generator = SyntheticBarcodedBAMGenerator(
        _test_num_cells,
        _test_max_genes,
        _test_annotation_file,
        _test_gene_expression_rate)

    with tempfile.TemporaryDirectory() as _test_temp_dir:
        # generate test data
        synthetic_data_generator.generate_synthetic_bam_and_counts_matrix(
            _test_temp_dir,
            _test_num_duplicates,
            _test_num_missing_some_tags,
            _test_num_multiple_gene_alignments,
            _test_max_gene_hits_per_multiple_gene_alignments,
            query_sort_order=query_sort_order)

        # test data paths
        test_bam_path = os.path.join(_test_temp_dir, SyntheticBarcodedBAMGenerator.bam_output_filename)
        test_count_matrix_path = os.path.join(_test_temp_dir, SyntheticBarcodedBAMGenerator.count_matrix_output_filename)
        test_row_index_path = os.path.join(_test_temp_dir, SyntheticBarcodedBAMGenerator.row_index_output_filename)
        test_col_index_path = os.path.join(_test_temp_dir, SyntheticBarcodedBAMGenerator.col_index_output_filename)

        # create CountMatrix from the synthetic bam
        count_matrix_from_bam: CountMatrix = CountMatrix.from_sorted_tagged_bam(test_bam_path, _test_annotation_file)

        # load the test counts matrix
        count_matrix_data_expected = np.load(test_count_matrix_path)
        row_index_expected = np.load(test_row_index_path)
        col_index_expected = np.load(test_col_index_path)

    count_matrix_data_from_bam = count_matrix_from_bam.matrix.todense()
    row_index_from_bam = count_matrix_from_bam.row_index
    col_index_from_bam = count_matrix_from_bam.col_index

    # sort expected and from_bam results by their respective row and column indices, since their sort order
    # is not part of the design specs and is considered arbitrary
    (sorted_count_matrix_data_from_bam,
     sorted_row_index_from_bam,
     sorted_col_index_from_bam) = _get_sorted_count_matrix(
        count_matrix_data_from_bam, row_index_from_bam, col_index_from_bam)
    (sorted_count_matrix_data_expected,
     sorted_row_index_expected,
     sorted_col_index_expected) = _get_sorted_count_matrix(
        count_matrix_data_expected, row_index_expected, col_index_expected)

    # assert equality
    assert np.allclose(sorted_count_matrix_data_from_bam, sorted_count_matrix_data_expected)
    assert all([l == r for l, r, in zip(sorted_row_index_from_bam, sorted_row_index_expected)])
    assert all([l == r for l, r, in zip(sorted_col_index_from_bam, sorted_col_index_expected)])
