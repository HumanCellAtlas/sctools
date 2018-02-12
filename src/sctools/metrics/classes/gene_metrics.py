from .base_metric import BaseMetric
from collections import Counter, defaultdict, namedtuple

# fragment could also just be a numpy structured array of length 5 with
# dtype (int32, int8, bool, bool, bool)
Fragment = namedtuple('Fragment', ['observations', 'min_mismatches', 'is_unique', 'is_spliced'])


class TagError(Exception):
    pass


class GeneMetrics(BaseMetric):
    """
    measure number of total genes

    For each gene:
    --------------
    number of reads
    number of fragments
    number of UMIs
    number of spliced and unspliced alignemnts
    number of wrong-strand alignments
    number of cells
    nubmer of substitution errors (mean)
    fraction of observations with single-read evidence
    fraction of observations
    fraction of reads mapped outside window

    higher order:
    -------------
    highest PC that the gene correlates with / loads onto
    highest DC that the gene correlates with / loads onto
    number of clusters gene is detected in

    these metrics more or less require a count matrix
    there are also requirements for a gtf file
    calculate the metrics that need neither first

    """

    def __init__(self, sample_name, cell_barcode_tag='CB', molecular_barcode_tag='UB'):
        self.basename = sample_name
        self.cell_barcode_tag = cell_barcode_tag
        self.molecular_barcode_tag = molecular_barcode_tag
        # dict[cell[gene[umi[fragment]]]] = Fragment()
        self._metric_data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

    def gather_metric(self, sam_record):
        """Extract metric information from a sam_record

        :param pysam.AlignedFragment sam_record: pysam object generated from a SAM format record
        """
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

        # record is aligned
        if sam_record.flag & 4:
            position, chromosome, strand, gene = 0, 0, 0, -2  # unaligned read
        else:
            # parse coding XF tag first, use this information to store the gene information
            try:
                alignment_type = sam_record.get_tag('XF')
                if alignment_type == 'CODING':
                    try:
                        gene = sam_record.get_tag('GE')
                    except KeyError:
                        raise TagError('gene not found with the expected tag: %s' % 'GE')
                elif alignment_type == 'INTRONIC':
                    gene = 0  # aligned, but to an intron
                else:
                    gene = -1  # aligned, but to integenic region
            except:
                raise TagError('should have a tag')

        # after here I kinda gave up and started pseudocoding
        # this is subtly wrong; needs to be done differently if not aligned.
        position = sam_record.position
        strand = not sam_record & 16  # 1 is +, 0 is -
        chromosome = sam_record.reference_id  # int instead of string, smaller
        min_mismatches = sam_record.cigar_tuples()  # do something with the cigar
        is_unique = sam_record.get_tag('')  # figure out the right tag
        is_spliced = sam_record.get_tag('')  # figure out the right tag

        try:
            this_fragment = self._metric_data[cell_barcode][gene][molecular_barcode][(position, chromosome, strand)]
        except KeyError:
            this_fragment = Fragment(1, min_mismatches, is_unique, is_spliced)

        raise NotImplementedError

    def calculate_and_output(self):
        raise NotImplementedError
