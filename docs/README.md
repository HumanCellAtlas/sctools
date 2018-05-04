# Build Docs 

1. Make sure you have [Sphinx](http://www.sphinx-doc.org/en/master/) installed.
2. Install the sctools package in advance following the instructions.
3. From the current directory, type:

```bash
make install
```

Note that there are still some bugs to be worked out. 
- There are warnings about: 
```
WARNING: [autosummary] failed to import 'sctools.metrics.CellMetrics': no module named sctools.metrics.CellMetrics
WARNING: [autosummary] failed to import 'sctools.metrics.GeneMetrics': no module named sctools.metrics.GeneMetrics
WARNING: [autosummary] failed to import 'sctools.metrics.MetricAggregatorBase': no module named sctools.metrics.MetricAggregatorBase
```

- There are a bunch of warnings: `WARNING: Unexpected section title.`
- There are a bunch of warnings: `WARNING: toctree contains reference to nonexisting document`

Most of the warnings can be solved by refactoring the docstrings and standardize the usages of `autosummary` later.