import abc


class BaseMetric:

    @abc.abstractmethod
    def setup(self):
        """Do some kind of one time initialization for the metrics class"""

    @abc.abstractmethod
    def accept_read(self, sam_record):
        """What the metric to do with each read of a bam"""

    @abc.abstractmethod
    def finish_up(self):
        """Clean up stuff / compute final metrics / write files"""
