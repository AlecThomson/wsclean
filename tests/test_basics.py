import pytest
import sys
from utils import validate_call

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf


@pytest.mark.parametrize("command", ["-version", "-this-is-not-a-valid-parameter"])
def test_basis(command):
    if command == "-this-is-not-a-valid-parameter":
        with pytest.raises(Exception):
            validate_call([tcf.WSCLEAN, command])
    else:
        validate_call([tcf.WSCLEAN, command])
