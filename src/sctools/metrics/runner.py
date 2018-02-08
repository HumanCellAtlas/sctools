import argparse
import pysam

class Runner:

    def __init__(self, sam_file, open_mode=None):
        """
        :param str alignment_file: sam or bam file.
        :param str open_mode: optional, mode to read file. Will be autodetected by file type if
          the file contains the correct suffix for its type.
        """

        if open_mode is None:
            if sam_file.endswith('.bam'):
                open_mode = 'rb'
            elif sam_file.endswith('.sam'):
                open_mode = 'r'
            else:
                raise ValueError('could not autodetect file type for alignment_file %s '
                                 '(detectible suffixes: .sam, .bam)' % sam_file)
        self.input_sam = sam_file
        self.open_mode = open_mode

    def run_metrics(self, metrics_to_run):
        """
        :param list metrics_to_run: list of metrics to run for each record in input_sam
        """

        for metric in metrics_to_run:
            metric.initialize()

        with pysam.AlignmentFile(self.input_sam, self.open_mode) as sam:
            for record in sam:
                for metric in metrics_to_run:
                    metric.gather_metric(record)

        for metric in metrics_to_run:
            metric.calculate_and_output()


def convert_class_name_to_class(metric_class_names):
    # maybe some kind of case statement or python equivalent to grab a class from a string
    # once we define a metric class that extends base_metric.py we can fill out this method
    pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_bam", help="path to input bam file")
    parser.add_argument("-m", "--metrics", nargs='+', help="metrics to run on input_bam")
    args = parser.parse_args()

    metric_classes = [collect_unique_fragment_by_umi.UniqueFragmentPerUMI()]  # convert_class_name_to_class(args.metrics)

    runner = Runner(args.input_bam)
    runner.run_metrics(metric_classes)


if __name__ == '__main__':
    main()
