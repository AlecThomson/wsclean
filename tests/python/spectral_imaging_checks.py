import pytest
import os
import sys
from utils import validate_call
from astropy.io import fits
import numpy

# Append current directory to system path in order to import testconfig
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


def name(name: str):
    return os.path.join(tcf.RESULTS_DIR, name)
