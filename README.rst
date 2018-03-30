sctools
=======

.. image:: https://circleci.com/gh/HumanCellAtlas/sctools/tree/master.svg?style=svg
    :target: https://circleci.com/gh/HumanCellAtlas/sctools/tree/master

sctools provides utilities for manipulating sequence data formats suitable for use in distributed
systems analyzing large biological datasets.

Download and Installation
-------------------------

.. code bash
   git clone https://github.com/humancellatlas/sctools.git
   cd sctools
   pip3 install .
   pytest  # verify installation; run tests

sctools Package
---------------

The sctools package provides both command line utilities and classes designed for use in python
programs.

Command Line Utilities
----------------------
1. Attach10XBarcodes: given an unaligned bam file, attach barcodes from an identically ordered
forward and forward index fastq file

Main Package Classes
--------------------

1. **Platform**: an abstract class that defines a common data structure for different 3' sequencing
   formats. All algorithms and methods in this package that are designed to work on 3' sequencing data
   speak to this common data structure. Currently 10X_v2 is defined.

2. **Reader**: a general iterator over arbitrarily zipped file(s) that is extended to work with common
   sequence formats like fastq (fastq.Reader) and gtf (gtf.Reader). We recommend using the pysam
   package for reading sam and bam files.

3. **TwoBit & ThreeBit** DNA encoders that store DNA in 2- and 3-bit form. 2-bit is smaller but
   randomizes "N" nucleotides. Both classes support fastq operations over common sequence tasks such
   as the calculation of GC content.

4. **ObservedBarcodeSet & PriorBarcodeSet**: classes for analysis and comparison of sets of barcodes
   such as the cell barcodes used by 10X genomics. Supports operations like summarizing hamming
   distances and comparing observed sequence diversity to expected (normally uniform) diversity.

5. **gtf.Reader & gtf.Record** GTF iterator and GTF record class that exposes the gtf
   fields as a lightweight, lazy-parsed python object.

6. **fastq.Reader & fastq.Record** fastq reader and fastq record class that exposes the fastq fields
   as a lightweight, lazy-parsed python object.


Viewing Test Results and Coverage
---------------------------------
to calculate and view test coverage cd to the ``sctools`` directory and
type the following two commands to generate the report and open it in your web browser:

.. code bash
   pytest --cov-report html:cov_html --cov=sctools
   open cov_html/index.html

Definitions
-----------

1. Fragment

2. Molecule

3. Cell Barcode

4. Molecule Barcode

5. Bam/Sam file: unless specified, refers to any aligned or unaligned bam/sam file
# todo add more
