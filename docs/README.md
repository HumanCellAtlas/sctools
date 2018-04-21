# Build Docs 
from the repository root, type: 

```bash
sphinx-apidoc -f -o docs/source/ ./
sphinx-apidoc -f -o docs/source/ src/sctools
sphinx-build -b html docs/source/ docs/build/
```

Note that there are still some bugs to be worked out. 
- `autosummary` is not a known directive
- `modules.rst` is not included in any toctree
- `setup` module cannot be imported
- `sctools.rst` is not included in any toctree
- /Users/ajc/projects/humancellatlas/sctools/src/sctools/reader.py:docstring of sctools.reader.zip_readers:: WARNING: more than one target found for cross-reference 'Reader': sctools.fastq.Reader, sctools.gtf.Reader, sctools.reader.Reader
- there are a bunch of unexpected section titles
