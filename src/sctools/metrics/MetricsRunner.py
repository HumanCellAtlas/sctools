import argparse


class MetricsRunner:

    @staticmethod
    def run_metrics(list_of_metrics_to_run, input_bam):

        for metric in list_of_metrics_to_run:
            metric.setup()

        for read in input_bam:
            metric.accept_read(read)

        for metric in list_of_metrics_to_run:
            metric.calculate_and_output()


def convert_class_name_to_class(list_of_class_names):
    # maybe some kind of case statement or python equivalent to grab a class from a string
    pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_bam", help="path to input bam file")
    parser.add_argument("-m", "--metric", action="append", help="metrics to run on input_bam")
    args = parser.parse_args()
    metric_classes = None  # convert_class_name_to_class(args.metric)
    MetricsRunner.run_metrics(metric_classes, args.input_bam)


if __name__ == '__main__':
    main()
