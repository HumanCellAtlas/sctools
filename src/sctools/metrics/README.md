# Loose structure of code

#### BaseMetric - base class for different metrics collectors
* setup() - any setup a metrics class has to do
* accept_read() - process a record for your metric
* finish_up() - finalize any calculations needed and write out file


#### MetricsRunner - class with command line interface, determines what metrics to run and inputs for those metrics
* loop through each metric and call setup() on each
* loop through bam calling metric.accept_read() for each metric
* loop through each metric and call finish_up() on each 
