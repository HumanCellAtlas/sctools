from . import BaseMetric


class CounterMetric(BaseMetric):
    count = 0

    def setup(self):
        print("do nothing")

    def accept_read(self, sam_record):
        global count
        if sam_record.isUnmapped():
            count += 1 # keep track of how many unmapped reads we've seen

    def finish_up(self):
        # output metrics for unmapped reads
