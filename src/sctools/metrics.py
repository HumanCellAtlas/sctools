from typing import Iterator
from collections import Counter


class GatherCellMetrics:

    # assumes file is sorted by cell barcode then molecule barcode

    def __init__(self):
        pass

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

    # need some kind of data structure to hold cell metrics
    def read_bam(self, bam_file: str):

        cell_barcode_histogram = Counter()
        # cell entropy
        # number cell barcodes
        # valid barcodes
        molecule_barcode_histogram = Counter()
        # umi entropy

        for cell in self.iter_cells():

            # maybe make all this stuff into a custom molecule with a .csv write function
            number_reads = 0
            number_molecules = 0
            # reads / molecule

            number_noise_reads = 0
            number_cell_barcode_bases_above_30 = 0
            number_perfect_cell_barcodes = 0  # inv: fraction cells with errors
            average_number_molecule_barcode_bases_above_30 = 0
            number_perfect_molecule_barcodes = 0  # inv: fraction molecules with errors
            average_number_genomic_reads_quality_above_30 = 0

            number_reads_mapped_to_UTR = 0
            number_reads_mapped_to_transcriptome = 0
            number_reads_mapped_intergenic = 0
            number_reads_mapped_intronic = 0
            number_reads_mapped_uniquely = 0
            number_reads_mapped_too_many_loci = 0
            number_reads_mapped_multiple = 0
            number_reads_mapped_outside_window = 0  # reads should be within 1000 bases of UTR
            number_duplicate_reads = 0
            number_genes_detected = 0

            plus_strand_reads = 0  # strand balance

            number_spliced_reads = 0
            number_antisense_reads = 0

            average_number_fragments_per_molecule = 0
            number_molecules_with_single_read_evidence = 0

            number_molecules_mapped_to_UTR = 0
            number_molecules_mapped_to_transcriptome = 0
            number_molecules_mapped_intergenic = 0
            number_molecules_mapped_intronic = 0

            # higher-order statistics:
            yield_ = None

            for molecule in self.iter_molecules(cell):
                pass  # fill in the data described above

        # post-process the histograms


    class GatherGeneMetrics():

        # assumes file is sorted by GE tag

        def __init__(self):
            pass

        def iter_genes(self, bam_file) -> Iterator:
            pass

        def read_bam(self, bam_file: str):

            for gene in self.iter_genes(bam_file):
                pass # fill in the data described above
