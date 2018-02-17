from .encodings import TwoBit
from collections import defaultdict
import argparse
import pysam
# from typing import NamedTuple


class TagError(Exception):
    pass


# class Molecule(NamedTuple):
#     min_mismatches: int
#     reverse_strand: int
#     is_spliced: bool  # logic: if any reads overlapping an intron are found, set to False
#     is_unique: bool = True
#     observations: int = 0


class Molecule:

    __slots__ = ['min_mismatches', 'reverse_strand', 'is_spliced', 'is_unique', 'observations']

    def __init__(self, min_mismatches, reverse_strand, is_spliced, is_unique=True, observations=0):
        self.min_mismatches = min_mismatches
        self.reverse_strand = reverse_strand
        self.is_spliced = is_spliced
        self.is_unique = is_unique
        self.observations = observations

    def __repr__(self):
        return 'obs: %d, spliced: %r, reverse: %r, mismatches: %d' % (
            self.observations, self.is_spliced, self.reverse_strand, self.min_mismatches
        )


class MetricCollector:

    def __init__(self, sample_name, cell_barcode_tag='CB', molecular_barcode_tag='UB'):
        self.basename = sample_name
        self.cell_barcode_tag = cell_barcode_tag
        self.molecular_barcode_tag = molecular_barcode_tag
        self._metric_data = defaultdict(dict)

    def gather_metric_data(self, sam_record):
        """Extract metric information from a sam_record

        :param pysam.AlignedFragment sam_record: pysam object generated from a SAM format record
        """

        # get barcodes
        try:
            cell_barcode = sam_record.get_tag(self.cell_barcode_tag)
        except KeyError:
            raise TagError('cell barcode not found with the expected tag %s' %
                           self.cell_barcode_tag)

        try:
            molecular_barcode = sam_record.get_tag(self.molecular_barcode_tag)
        except KeyError:
            raise TagError('molecular barcode not found with the expected tag %s' %
                           self.molecular_barcode_tag)

        # convert barcodes to ints to reduce storage footprint
        # todo change encoder to take strings too
        cell_encoded = TwoBit.encode(cell_barcode.encode())
        molecule_encoded = TwoBit.encode(molecular_barcode.encode())
        try:
            is_multimapped = sam_record.get_tag('NH') > 1
        except KeyError:
            raise TagError('NH tag not found. This tag indicates the number of times a read '
                           'aligned. Is this bam file aligned?')

        # if the record is not mapped, increment and return
        if sam_record.is_unmapped:
            try:
                self._metric_data[(cell_encoded, molecule_encoded)]['umapped'] += 1
            except KeyError:
                self._metric_data[(cell_encoded, molecule_encoded)]['umapped'] = 1

        # if the record is multimapped, increment and return # todo think more about this
        elif is_multimapped:
            try:
                self._metric_data[(cell_encoded, molecule_encoded)]['multi-mapped'] += 1
            except KeyError:
                self._metric_data[(cell_encoded, molecule_encoded)]['multi-mapped'] = 1

        # record is uniquely mapped.
        else:
            # assumptions:
            # assume XF tag will contain all necessary information about where a gene is mapped
            # open question about how to store reads that overlap splice junctions that will
            # need to be treated in another module (TagGeneExon), here, assume that unspliced reads
            # that contain introns are labeled as intronic.
            # assume such genes will also have a GE tag (this may not be the case, in which case
            # we should correct it)
            try:
                alignment_type = sam_record.get_tag('XF')
            except KeyError:
                raise TagError('No XF tag found. All aligned records should have an XF tag '
                               'declaring the part of the genome they are aligned against. Options '
                               'are {CODING, INTRONIC, INTERGENIC}')
            if alignment_type == 'INTERGENIC':
                try:
                    self._metric_data[(cell_encoded, molecule_encoded)]['intergenic'] += 1
                except KeyError:
                    self._metric_data[(cell_encoded, molecule_encoded)]['intergenic'] = 1
            else:
                try:
                    gene = sam_record.get_tag('GE')
                except KeyError:
                    try:
                        self._metric_data[(cell_encoded, molecule_encoded)]['no-gene'] += 1
                    except KeyError:
                        self._metric_data[(cell_encoded, molecule_encoded)]['no-gene'] = 1
                    return  # TagGeneExon not adding gene tags where it should be
                    # raise TagError('No GE tag found. Any read aligned to CODING or INTERGENIC '
                    #                'sequence should have a GE tag indicating the gene body it '
                    #                'aligned against. Read:\n%s' % str(sam_record))
                n_nucleotides, n_blocks = sam_record.get_cigar_stats()
                min_mismatches = sam_record.infer_query_length() - n_nucleotides[0]
                is_unique = sam_record.get_tag('NH') == 1
                is_spliced = n_nucleotides[3] > 0  # todo should this number be larger?
                # want to know if read is outside window (stretch, add later)
                try:
                    molecule = self._metric_data[(cell_encoded, molecule_encoded)][gene]
                    molecule.observations += 1
                    molecule.is_unique &= is_unique
                    molecule.is_spliced &= is_spliced
                    molecule.min_mismatches = min(molecule.min_mismatches, min_mismatches)
                    if sam_record.is_reverse:
                        molecule.reverse_strand += 1
                except KeyError:
                    reverse_strand = 1 if sam_record.is_reverse else 0
                    self._metric_data[(cell_encoded, molecule_encoded)][gene] = Molecule(
                        min_mismatches=min_mismatches, is_unique=is_unique, is_spliced=is_spliced,
                        observations=1, reverse_strand=reverse_strand
                    )

    # metric calculators
    def calculate_and_output(self):
        """Clean up stuff / compute final metrics / write files"""
        raise NotImplementedError

    def _noop(self):
        pass  # internal, for testing


class Runner:

    def __init__(self, sam_file, cell_barcode_tag, molecule_barcode_tag, basename, open_mode=None):
        """
        :param str sam_file: sam or bam file.
        :param str open_mode: optional, mode to read file. Will be autodetected by file type if
          the file contains the correct suffix for its type.
        """

        if open_mode is None:
            if sam_file.endswith('.bam'):
                open_mode = 'rb'
            elif sam_file.endswith('.sam'):
                open_mode = 'r'
            else:
                raise ValueError('could not autodetect file type for sam_file %s '
                                 '(detectible suffixes: .sam, .bam)' % sam_file)
        self.input_sam = sam_file
        self.open_mode = open_mode
        self.metric_accumulator = MetricCollector(
            basename, cell_barcode_tag, molecule_barcode_tag)

    def run_metrics(self, metrics_to_run):

        if not isinstance(metrics_to_run, (list, tuple)):
            raise TypeError('metrics must be a list of tuple of metrics')

        with pysam.AlignmentFile(self.input_sam, self.open_mode) as sam:
            for record in sam:
                self.metric_accumulator.gather_metric_data(record)

        for metric in metrics_to_run:
            metric_name = getattr(self.metric_accumulator, metric)
            metric_name()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_bam", help="path to input bam file")
    parser.add_argument("basename", type=str, help="basename of metrics to output")
    # todo can switch this to "choices"
    parser.add_argument(
        "-m", "--metrics", nargs='+',
        help='metrics to run on input_bam, possible values include: UniqueFragmentPerUMI.')
    parser.add_argument(
        "-c", "--cell-barcode-tag", type=str,
        help="tag value to grab cell barcode from record ", default="CR")
    parser.add_argument(
        "-m", "--molecular-barcode-tag", type=str,
        help="tag value to grab cell barcode from record", default="UR")
    args = parser.parse_args()

    if len(args.metrics) == 0:
        print("you need to provide at least one metric for this program to run")
        exit(1)

    runner = Runner(args.input_bam, args.metrics)
    runner.run_metrics(args.cell_barcode_tag, args.molecular_barcode_tag, args.basename,
                       args.metrics)


if __name__ == '__main__':
    main()
