"""
Global constants
================

.. currentmodule:: sctools

This module contains global constants, such as various barcoded BAM tags, and sctools-specific
constants.
"""

# BAM tag constants

RAW_SAMPLE_BARCODE_TAG_KEY = 'SR'
QUALITY_SAMPLE_BARCODE_TAG_KEY = 'SY'

MOLECULE_BARCODE_TAG_KEY = 'UB'
RAW_MOLECULE_BARCODE_TAG_KEY = 'UR'
QUALITY_MOLECULE_BARCODE_TAG_KEY = 'UY'

CELL_BARCODE_TAG_KEY = 'CB'
RAW_CELL_BARCODE_TAG_KEY = 'CR'
QUALITY_CELL_BARCODE_TAG_KEY = 'CY'

GENE_NAME_TAG_KEY = 'GE'
NUMBER_OF_HITS_TAG_KEY = 'NH'

ALIGNMENT_LOCATION_TAG_KEY = 'XF'
INTRONIC_ALIGNMENT_LOCATION_TAG_VALUE = 'INTRONIC'
CODING_ALIGNMENT_LOCATION_TAG_VALUE = 'CODING'
UTR_ALIGNMENT_LOCATION_TAG_VALUE = 'UTR'
INTERGENIC_ALIGNMENT_LOCATION_TAG_VALUE = 'INTERGENIC'

INDROP_SPLITTER_SEQUENCE = 'GAGTGATTGCTTGTGACGCCTT'

def generate_splitter_regex() :
    splitters = []
    search_substring = INDROP_SPLITTER_SEQUENCE[:6]
    characters = ['A', 'T', 'C', 'G']
    for x in range(0, 6):
        splitters.extend(search_substring[:x] + c + search_substring[x:] for c in characters)
    regex_string = '(' + '|'.join(splitters) + ')'
    return regex_string


SPLITTER_REGEX = generate_splitter_regex()

# search for splitters

# bam.py constants

MAX_BAM_SPLIT_SUBFILES_TO_WARN = 500
MAX_BAM_SPLIT_SUBFILES_TO_RAISE = 1000
