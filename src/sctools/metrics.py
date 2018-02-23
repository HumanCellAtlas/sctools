from typing import Iterator
from collections import Counter


class CellBarcodeMetrics:

    def __init__(self):

        # count information
        self._number_reads = 0
        self._number_molecules = 0
        # dependent metric: reads / molecule

        # barcode quality data
        self._number_noise_reads = 0  # long polymers, N-sequences
        self._number_cell_barcode_bases_above_30 = 0
        self._number_perfect_cell_barcodes = 0  # inv: fraction cells with errors
        self._average_number_molecule_barcode_bases_above_30 = 0
        self._number_perfect_molecule_barcodes = 0  # inv: fraction molecules with errors
        self._average_number_genomic_reads_quality_above_30 = 0

        # alignment locations
        self._number_reads_mapped_to_UTR = 0
        self._number_reads_mapped_to_transcriptome = 0
        self._number_reads_mapped_intergenic = 0
        self._number_reads_mapped_intronic = 0
        self._number_reads_mapped_uniquely = 0
        self._number_reads_mapped_too_many_loci = 0
        self._number_reads_mapped_multiple = 0
        self._number_reads_mapped_outside_window = 0  # reads should be within 1000 bases of UTR
        self._number_duplicate_reads = 0
        self._number_genes_detected = 0

        self._plus_strand_reads = 0  # strand balance

        self._number_spliced_reads = 0
        self._number_antisense_reads = 0

        self._average_number_fragments_per_molecule = 0
        self._number_molecules_with_single_read_evidence = 0

        self._number_molecules_mapped_to_UTR = 0
        self._number_molecules_mapped_to_transcriptome = 0
        self._number_molecules_mapped_intergenic = 0
        self._number_molecules_mapped_intronic = 0

        # higher-order statistics:
        self._yield = 0


class GatherCellBarcodeMetrics:

    def __init__(self, bam_file: str):
        self._bam_file = bam_file

    def iter_cells(self) -> Iterator:
        """
        function to iterate over all the cells of a bam file
        """
        raise NotImplementedError

    def iter_molecules(self, cell) -> Iterator:
        """
        function to iterate over all the molecules of a cell
        """
        raise NotImplementedError

    def read_bam(self):

        cell_barcode_histogram = Counter()
        # cell entropy
        # number cell barcodes
        # valid barcodes

        molecule_barcode_histogram = Counter()
        # umi entropy

        for cell in self.iter_cells(self._bam_file):

            # maybe make all this stuff into a custom cell class with a .csv write function

            for molecule in self.iter_molecules(cell):
                pass  # fill in the data described above

        # post-process the histograms for the metrics they contain


class Gene():

    def __init__(self):

        # read metrics
        self._number_reads_mapped_to_UTR = 0
        self._number_reads_mapped_to_transcriptome = 0
        self._number_reads_mapped_intronic = 0
        self._number_reads_mapped_uniquely = 0
        self._number_reads_mapped_multiple = 0
        self._number_reads_mapped_outside_window = 0  # reads should be within 1000 bases of UTR
        self._number_duplicate_reads = 0
        self._average_distance_from_termination_site = 0

        # strand and splicing
        self._plus_strand_reads = 0  # strand balance
        self._number_spliced_reads = 0
        self._number_antisense_reads = 0

        # molecule confidence measurements
        self._average_number_fragments_per_molecule = 0
        self._number_molecules_with_single_read_evidence = 0

        # molecule metrics
        self._number_molecules_mapped_to_UTR = 0
        self._number_molecules_mapped_to_transcriptome = 0
        self._number_molecules_mapped_intronic = 0

        # higher-order statistics:
        self._yield = 0
        self._exon_histogram = Counter()


class GatherGeneMetrics:

    # assumes file is sorted by GE tag

    def __init__(self):
        pass

    def iter_genes(self, bam_file) -> Iterator:
        raise NotImplementedError

    def read_bam(self, bam_file: str):

        for gene in self.iter_genes(bam_file):
            pass  # fill in the data described above

