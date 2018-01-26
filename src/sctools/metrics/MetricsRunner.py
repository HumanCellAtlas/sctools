class MetricsRunner:
    list_of_base_metrics_to_run = []
    input_bam = "some/path/through/command"
    """
    Based on input passed through the command line we would decide
    which metrics to run on what bam
    """

    def run_metrics(self):
        global list_of_base_metrics_to_run, input_bam

        for metric in list_of_base_metrics_to_run:
            metric.setup()

        for read in input_bam:
            metric.accept_read(read)

        for metric in list_of_base_metrics_to_run:
            metric.finish_up()
