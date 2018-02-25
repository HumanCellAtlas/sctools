## Overall Metric Processing
Metric processing relies on the structure of the Optimus pipeline. Specifically, it:

1. relies upon chunking the input bam files into smaller pieces
2. relies on the guarantee of the optimus pipeline that for each cell, all reads associated with it
are found in the same chunk
3. Relies on sorting the (small) bam files based on different tags by first converting them to sam. 
The distributed nature of the pipeline makes this practical.
4. This sorting enables us to process all of the data for a molecule at one time. This is important, 
as many metrics integrate over a molecule's data. Without sorting, a data structure must be built to 
build up the information from the bam file, and the resulting memory footprint is 20-40gb for a 
standard HiSeq run.  
5. To sort, the file needs three tags: Cell barcode (CB), Molecule barcode (UB) and GeneMetrics (GE). The 
latter tag is not, to my knowledge, in the HTSlib spec. 

Despite the reliance on optimus assumptions, We will be able to release this tool as a stand alone 
that can run on an arbitrary, properly tagged bam file, but it will rely on other sctools 
functionality and may NOT be effectively portable across operating systems (relies on linux sort)

An abbreviated WDL could be made in the future which would contain: 

```
1. SplitBamByCellBarcode
2. scatter[CalculateMetrics]
3. MergeMetrics
```

## Implementation: 

I believe the best course of action is to generate a metrics file that aligns with our eventual 
count matrix. This is consistent with scanpy/loom formats which use DataFrame-like objects to store
metadata for both samples and features. It will enable the most common access patterns of comparing
metrics across the relevant data, and more complex access patterns that fix mixed-effect models to 
correct for metadata values. 

This means we need two metrics outputs: one for Cells, and one for genes. 

I've pseudo-coded three classes, each of which must be subclassed to write a Metric: 
```
SequenceMetricAggregator:
  - CellMetricAggregator
  - GeneMetricAggregator

MetricGatherer:
  - CellMetricGatherer
  - GeneMetricGatherer
```

SequenceMetricAggregator stores the information for a unit of the relevant data (cell, gene) as records are
being processed, and has csv writing capability so that the information for only one unit is stored
at a time, ensuring scalability to extremely large datasets. It has a `parse_records()` function 
which has not yet been implemented which will do most of the work

MetricGatherer defines iterables to do the record grouping, the bam parsing pattern necessary, 
and a `check_sort_order` function to verify the file is in the required order. 

MergeMetrics merges multiple metric outputs from the scattered chunks. This is a trivial 
concatenation in the case of cell metrics, and a more complex merge in the case of gene metrics. The
`merge_metrics` function is not yet implemented, but will carry this out. 


