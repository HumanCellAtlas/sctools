from . import bam
from . import encodings
from . import barcode
from . import fastq
from . import gtf
from . import stats
from . import reader
from . import metrics
from . import platform
from . import consts
from . import groups
from pkg_resources import get_distribution, DistributionNotFound


try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    pass
