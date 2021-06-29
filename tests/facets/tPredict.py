import os
import pytest
from subprocess import check_call, check_output
import shutil
import logging

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf

# Prepend path with current working directory to make sure
# wsclean executable from the build directory
os.environ["PATH"] = f"{os.getcwd()}:{os.environ['PATH']}"

MWA_MOCK_ARCHIVE="MWA_ARCHIVE.tar.bz2"
MWA_MOCK_MS="MWA_MOCK.ms"
MWA_MOCK_FULL="MWA_MOCK_FULL.ms"
MWA_MOCK_FACET="MWA_MOCK_FACET.ms"


@pytest.fixture(autouse=True)
def prepare_ms():
    # Change to directory containing the data
    os.makedirs(tcf.WORKDIR, exist_ok=True)
    os.chdir(tcf.WORKDIR)
    if not os.path.isfile(MWA_MOCK_ARCHIVE):
        wget = f"wget -q www.astron.nl/citt/EveryBeam/MWA-single-timeslot.tar.bz2 -O {MWA_MOCK_ARCHIVE}"
        check_call(wget.split())

    os.makedirs(MWA_MOCK_MS, exist_ok=True)
    check_call(f"tar -xf {MWA_MOCK_ARCHIVE}  -C {MWA_MOCK_MS} --strip-components=1".split())

    # Not pythonic, but works fine and fast
    check_call(f"cp -r {MWA_MOCK_MS} {MWA_MOCK_FULL}".split())
    check_call(f"cp -r {MWA_MOCK_MS} {MWA_MOCK_FACET}".split())

@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    # Fixture is run only at end of test session
    # NOTE: to enable clean-up of the
    # MWA ms and the coefficients file, make sure
    # CLEANUP is in your environment variables

    # Remove measurement set
    def remove_mwa():
        if "CLEANUP" in os.environ:
            os.chdir(tcf.WORKDIR)
            os.remove(MWA_MOCK_ARCHIVE)
            shutil.rmtree(os.environ["MWA_MOCK_FULL"])
            shutil.rmtree(os.environ["MWA_MOCK_FACET"])

    request.addfinalizer(remove_mwa)


def assert_taql(command, expected_rows=0):
    taqlexe = shutil.which("taql")
    if taqlexe is not None:
        result = check_output([taqlexe, "-noph", command]).decode().strip()
        assert result == f"select result of {expected_rows} rows", result
    else:
        raise OSError("taql executable not found. Make sure the taql executable - provided by casacore(-tools) - is on your path.")


def predict_full_image(ms, gridder):
    # Predict full image
    s = f"wsclean -predict {gridder} -name point-source {ms}"
    check_call(s.split())

def predict_facet_image(ms, gridder):
    # Predict full image
    s = f"wsclean -predict {gridder} -facet-regions {tcf.FACETFILE_4FACETS} -name point-source {ms}"
    check_call(s.split())

# FIXME: we should test wstacking and wgridder here too
# but the fail on the taql assertion
@pytest.mark.parametrize('gridder', ["-use-wgridder"])
def test_predict(gridder):
    predict_full_image(MWA_MOCK_FULL, gridder)
    predict_facet_image(MWA_MOCK_FACET, gridder)
    taql_command = f"select from {MWA_MOCK_FULL} t1, {MWA_MOCK_FACET} t2 where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,5e-3))"
    assert_taql(taql_command)
