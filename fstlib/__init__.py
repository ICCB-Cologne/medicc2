import os
import sys
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

IS_WINDOWS = os.name=="nt"

class Semiring:
    REAL = 'real'
    LOG = 'log'
    TROPICAL = 'standard'

DEF_GAP_SYMBOL = '-'
DEF_DELTA = 0.0009765625
MAX_INT32 = 2**31 - 1

from fstlib.cext.pywrapfst import *
from fstlib.core import *
from fstlib.ext import *
import fstlib.algos
import fstlib.factory
import fstlib.tools




