import pytest
import pytest_lazyfixture
# from subprocess import check_call
import os

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf

# MWA_MOCK_ARCHIVE = "MWA_ARCHIVE.tar.bz2"
# MWA_MOCK_MS = "MWA_MOCK.ms"

# @pytest.fixture
def gridders():
    return {"wstakcing": "", "wgridder": "-use-wgridder", "idg": "-use-idg"}

# @pytest.fixture(scope="session", autouse=True)
# def prepare():
#     os.makedirs(tcf.WORKDIR, exist_ok=True)
#     os.chdir(tcf.WORKDIR)

#     if not os.path.isfile(MWA_MOCK_ARCHIVE):
#         wget = f"wget -q www.astron.nl/citt/EveryBeam/MWA-single-timeslot.tar.bz2 -O {MWA_MOCK_ARCHIVE}"
#         check_call(wget.split())

#     os.makedirs(MWA_MOCK_MS, exist_ok=True)
#     check_call(
#         f"tar -xf {MWA_MOCK_ARCHIVE}  -C {MWA_MOCK_MS} --strip-components=1".split()
#     )

#     os.path.isfile(tcf.FACETFILE_2FACETS)



# Test assumes that IDG and EveryBeam are installed
@pytest.mark.parametrize("gridder", gridders().items())
def test_stitching(gridder):
    print(gridder)
    print(type(gridder))
    # s = "wsclean "
    prefix = f"facet-stitch-{gridder[0]}"
    s = f"wsclean -quiet -size 256 256 -scale 4amin -pol XX,YY -facet-regions {tcf.FACETFILE_2FACETS} -name",  {prefix},  kMWA_MS


    pass
