import argparse
import pysam
import classes


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

    def run_metrics(self, args, metrics_to_run):
        """
        :param list metrics_to_run: list of metrics to run for each record in input_sam
        :param object args: all arguments passed to command line in order to pass to each metric class
        """

        for metric in metrics_to_run:
            metric.initialize(metric, args)

        with pysam.AlignmentFile(self.input_sam, self.open_mode) as sam:
            for record in sam:
                for metric in metrics_to_run:
                    metric.gather_metric(metric, record)

        for metric in metrics_to_run:
            metric.calculate_and_output(metric)


def convert_class_name_to_class(metric_class_names):
    return [getattr(classes, metric_name) for metric_name in metric_class_names]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_bam", help="path to input bam file")
    parser.add_argument("basename", type=str, help="basename of metrics to output")
    parser.add_argument("-m", "--metrics", nargs='+', help="metrics to run on input_bam")
    parser.add_argument("-cbt", "--cell_barcode_tag", type=str, help="tag value to grab cell barcode from record",
                        default="CR")
    parser.add_argument("-mbt", "--molecular_barcode_tag", type=str, help="tag value to grab cell barcode from record",
                        default="UR")
    args = parser.parse_args()

    if len(args.metrics) == 0:
        print("you need to provide at least one metric for this program to run")
        exit(1)

    metric_classes = convert_class_name_to_class(args.metrics)

    runner = Runner(args.input_bam)
    runner.run_metrics(args, metric_classes)


if __name__ == '__main__':
    main()
