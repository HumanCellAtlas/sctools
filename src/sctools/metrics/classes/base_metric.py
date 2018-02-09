class BaseMetric:

    def initialize(self, args):
        """Do some kind of one time initialization for the metrics class"""
        raise NotImplementedError

    def gather_metric(self, sam_record):
        """What the metric to do with each read of a bam
        :param Object sam_record: Object containing information from one read in a Sam file
        """
        raise NotImplementedError

    def calculate_and_output(self):
        """Clean up stuff / compute final metrics / write files"""
        raise NotImplementedError
