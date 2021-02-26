import os
import sys
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

IS_WINDOWS = os.name=="nt"

if __package__ is None:
    MODULE_BASE_PATH = "fstlib"
else:
    MODULE_BASE_PATH = os.path.dirname(sys.modules[__package__].__file__)

class Semiring:
    REAL = 'real'
    LOG = 'log'
    TROPICAL = 'standard'

DEF_GAP_SYMBOL = '-'
DEF_DELTA = 0.0009765625
MAX_INT32 = 2**31 - 1

try: 
    from fstlib.cext.pywrapfst import *
    log.info('Using internal openfst wrapper.')
except ImportError:
    try:
        from openfst_python import *
        log.info('Using openfst-python.')
    except ImportError:
        from pywrapfst import *
        log.info('Using openfst supplied pywrapfst.')

import fstlib
from fstlib.core import *
from fstlib.ext import *
import fstlib.algos
import fstlib.factory
import fstlib.tools




