import os

try:
    MS = os.environ["MS"]
except:
    raise IOError("Make sure MS environment variable is set, specifying the path to the measuremnt set that you want to use in the tests")

try:
    MWA_PATH = os.environ["MWA_PATH"]
except:
    MWA_PATH = ""


DIMS = "-size 1024 1024 -scale 1amin"
RECTDIMS = "-size 1536 1024 -scale 1amin"
FACETFILE_2FACETS = "../tests/data/ds9_2facets.reg"
FACETFILE_NFACETS = "../tests/data/ds9_nfacets.reg"