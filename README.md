#### MetricCollector - base class for different metrics collectors
* initialize() - any setup a metrics class has to do
* gather_metric() - process a record for your metric
* calculate_and_output() - finalize any calculations needed and write out file


#### MetricsRunner - class with command line interface, determines what metrics to run and inputs for those metrics
* loop through each metric and call initialize() on each
* loop through bam calling metric.gather_metric() for each metric
* loop through each metric and call calculate_and_output() on each

##### Command Line Usage
usage: runner.py [-h] [-m METRICS [METRICS ...]] [-cbt CELL_BARCODE_TAG]
                 [-mbt MOLECULAR_BARCODE_TAG]
                 input_bam basename
                 
Required arguments: input_bam, basename, and at least one metric


#### what kinds of objects do we need to parse
- sequence data
- quality data
- categorical tags
- cigar string
- bit flag

#### what do we need to do to the objects?
##### sequence data:
- count bases at each position
- count unique sequences (reads per molecule)
- count bases in sequence objects
- count errors in bases (cigar)
- calculate entropy (full sequence) -- might need to do this on keys (afterwards) or on sequences (during gather)
- partition reads with and without errors
- partition reads with and without multiple alignments

##### quality data:
- quality histogram by base
- quality histogram by nucleotide

##### categorical data:
- count number (could include genes)

#### what groupings/gather methods do we need?
- cell [cell]
- molecules [cell, umi, sample]
- sequence fragments [cell, umi, sample, [chrom, insert position] | sequence]
- gene, exon, transcript etc [_ tag]  # varies, but is just a single tag type, should be straightforward.
- alignments [qname | tag] might want to summarize multi-alignments in some way.
