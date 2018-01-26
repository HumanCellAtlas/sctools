from . import BaseMetric
import collections


class CounterMetric(BaseMetric):

    def setup(self):
        global c
        c = collections.Counter()

    def accept_read(self, sam_record):
        blah = sam_record.get_tag("CB") # grab tag value from input to main program and pass to here
        c[blah] += 1

    def finish_up(self):
        # output metric to a file
