from . import base_metric
import csv
import collections
import statistics


class UniqueFragmentPerUMI(base_metric.BaseMetric):

    def __init__(self):
        self.cell_barcode_to_umi = None
        self.basename = None
        self.cell_barcode_tag = None
        self.umi_barcode_tag = None

    def initialize(self, args):
        self.cell_barcode_to_umi = collections.defaultdict(dict)
        self.basename = args.basename
        self.cell_barcode_tag = args.cell_barcode_tag
        self.umi_barcode_tag = args.molecular_barcode_tag
        return

    def gather_metric(self, sam_record):
        if sam_record.is_unmapped:
            # there is a good possibility that we would want to track these metrics for unmapped reads as well
            # but as a first pass we will only do aligned reds
            return
        try:
            cell_barcode = sam_record.get_tag(self.cell_barcode_tag)
            umi_barcode = sam_record.get_tag(self.umi_barcode_tag)
        except KeyError:
            print("should we throw a custom error here or print out a warning/error to a log")
            return
        if sam_record.is_duplicate:
            return
        if umi_barcode in self.cell_barcode_to_umi[cell_barcode]:
            self.cell_barcode_to_umi[cell_barcode][umi_barcode] += 1
        else:
            self.cell_barcode_to_umi[cell_barcode][umi_barcode] = 1
        return

    def calculate_and_output(self):
        detail_header_columns = ["Cell_Barcode", "Total_UMI_Fragments", "Mean_Fragments_Per_UMI", "Standard_Deviation"]
        summary_header_columns = ["Total_UMI_Fragments", "Mean_Fragments_Per_UMI", "Standard_Deviation"]
        with open(self.basename + '.unique_fragment_by_umi_detail_metrics.tsv', 'w') as detail_metrics:
            writer = csv.writer(detail_metrics, delimiter='\t')
            writer.writerow(detail_header_columns)
            total_umi_fragments_captured_in_bam = 0
            total_umis_in_bam = 0
            fragments_per_umi_in_bam = []
            for cell_barcode in self.cell_barcode_to_umi.keys():
                total_umi_fragments_captured_in_cell_barcode = 0
                total_umis_in_cell_barcode = 0
                fragments_per_umi_in_cell_barcode = []
                for umi in self.cell_barcode_to_umi[cell_barcode].keys():
                    total_umis_in_cell_barcode += 1
                    fragments_for_umi = self.cell_barcode_to_umi[cell_barcode][umi]
                    total_umi_fragments_captured_in_cell_barcode += fragments_for_umi
                    fragments_per_umi_in_cell_barcode.append(fragments_for_umi)
                std_dev = 0.0
                if len(fragments_per_umi_in_cell_barcode) > 1:
                    std_dev = statistics.stdev(fragments_per_umi_in_cell_barcode)
                total_umi_fragments_captured_in_bam += total_umi_fragments_captured_in_cell_barcode
                total_umis_in_bam += total_umis_in_cell_barcode
                fragments_per_umi_in_bam.extend(fragments_per_umi_in_cell_barcode)
                writer.writerow([cell_barcode, total_umi_fragments_captured_in_cell_barcode,
                                 total_umi_fragments_captured_in_cell_barcode / total_umis_in_cell_barcode, std_dev])
        with open(self.basename + '.unique_fragment_by_umi_summary_metrics.tsv', 'w') as summary_metrics:
            writer = csv.writer(summary_metrics, delimiter='\t')
            writer.writerow(summary_header_columns)
            std_dev = statistics.stdev(fragments_per_umi_in_bam)
            writer.writerow([total_umi_fragments_captured_in_bam,
                             total_umi_fragments_captured_in_bam / total_umis_in_bam, std_dev])
        return

