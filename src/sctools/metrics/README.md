## Metric Processing
This module implements a metric suite that generates information on data quality at the level of 
both cells and genes. This QC information aligns with the cells and genes that make up the 
expression matrix, providing easy access to information that the user can examine to make decisions
about which cells or genes are of adequate quality to include in downstream processing. 

Metric processing in sctools can be run on large individual files, but also implements a map-reduce 
architecture execution at production scale. Specifically, the workflow is as follows: 

1. Chunk the input bam file using `SplitBam`, which generates several chunks, each of which is
guaranteed to contain all data for any cell it contains
2. Sort each chunk by cell, gene, and molecule tags to ensure that all the reads associated with 
a molecule are stored sequentially by cell (`CalculateCellMetrics`) or by gene 
(`CalculateGeneMetrics`)
3. For each cell or gene, parse the information by molecule, which typically loads fewer than 
10,000 records into memory at a time. 
4. Merge data across chunks using `MergeCellMetrics` or `MergeGeneMetrics`.

This map-reduce approach is currently implemented by the 
[HCA 3' pipeline](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl), 
but an abbreviated WDL could be made in the future which would contain: 

```
1. SplitBamByCellBarcode
2. scatter[CalculateMetrics]
3. MergeMetrics
```

## Implementation Details: 

This module implements 4 base classes that carry out metric processing. These are: 

```
SequenceMetricAggregator:
  - CellMetricAggregator
  - GeneMetricAggregator

MetricGatherer:
  - CellMetricGatherer
  - GeneMetricGatherer
 
MetricCSVWriter

MergeMetrics:
  - MergeCellMetrics
  - MergeGeneMetrics
```
MetricGatherer defines generator functions to group records into molecules, the bam parsing pattern 
necessary to process data iteratively. 

SequenceMetricAggregator stores the information for a unit of the relevant data (cell, gene), 
and processses all the records with the `.parse_records()` method. 

When all records of a single unit (cell, gene) have been processed, `.finalize()` is called to 
calculate any higher-order metrics (for example, the variance in quality scores across reads of the 
cell or gene), and it is written to file by `MetricSCVWriter`.  

MergeMetrics merges multiple metric outputs from the scattered chunks. This is a trivial 
concatenation in the case of cell metrics, and a more complex merge in the case of gene metrics. 
