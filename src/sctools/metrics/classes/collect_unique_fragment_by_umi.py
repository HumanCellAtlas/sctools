from .. import base_metric
import pysam


class UniqueFragmentPerUMI(base_metric):

    def __init__(self):
        self.cell_barcode_to_umi = {}
        self.umi_to_fragment = {}

    def initialize(self):
        return

    def gather_metric(self, sam_record):
        print("hi")
        return

    def calculate_and_output(self):
        print("done")
        return

