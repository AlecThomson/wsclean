import pytest
import sys
from utils import validate_call

# Append current directory to system path in order to import testconfig
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


@pytest.mark.parametrize("command", ["-version", "-this-is-not-a-valid-parameter"])
def test_basis(command):
    if command == "-this-is-not-a-valid-parameter":
        with pytest.raises(Exception):
            validate_call([tcf.WSCLEAN, command])
    else:
        validate_call([tcf.WSCLEAN, command])


@pytest.mark.parametrize("dd_psf_args", ["-dd-psf-grid 5 5", ""])
def test_dd_psfs_call(dd_psf_args):
    s = (
        f"{tcf.WSCLEAN} -size 200 200 -scale 2arcsec -niter 1 -mgain 0.7 -log-time -nmiter 1 -threshold 0.001 -use-idg "
        + dd_psf_args
        + f" {tcf.MWA_MOCK_MS}"
    )
    validate_call(s.split())
