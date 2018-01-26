from . import BaseMetric


class CounterMetric(BaseMetric):
    count = 0

    def setup(self):
        print("do nothing")

    def accept_read(self, sam_record):
        global count
        if sam_record.isUnmapped():
            count += 1

    def finish_up(self):
        # output histogram
