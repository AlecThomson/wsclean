import pytest
import shutil
import sys
import os
from subprocess import check_call

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf

# MODEL_IMAGE = "point-source-model.fits"
# MWA_COEFF_ARCHIVE = "mwa_full_embedded_element_pattern.tar.bz2"
# EVERYBEAM_BASE_URL = "http://www.astron.nl/citt/EveryBeam/"
# MWA_MOCK_MS = "MWA_MOCK.ms"
# MWA_MOCK_FULL = "MWA_MOCK_FULL.ms"
# MWA_MOCK_FACET = "MWA_MOCK_FACET.ms"


@pytest.fixture(scope="session", autouse=True)
def prepare():
    # Change to directory containing the data
    os.makedirs(tcf.WORKDIR, exist_ok=True)
    os.chdir(tcf.WORKDIR)
    os.makedirs(tcf.RESULTS_DIR, exist_ok=True)

    # Download and untar beam pattern file
    if not os.path.isfile(tcf.MWA_COEFF_FILE):
        check_call(["wget", "-q", tcf.EVERYBEAM_BASE_URL + tcf.MWA_COEFF_ARCHIVE])
        check_call(["tar", "xf", tcf.MWA_COEFF_ARCHIVE])


@pytest.fixture(scope="class")
def prepare_integration_tests():
    if not os.path.isfile(f"{tcf.MWA_MS}/table.f1"):
        check_call(["wget", "-q", os.path.join(EVERYBEAM_BASE_URL, f"{MWA_MS}.tgz")])
        check_call(["tar", "-xf", f"{tcf.MWA_MS}.tgz"])
        os.remove(f"{tcf.MWA_MS}.tgz")
    else:
        pass


@pytest.fixture(scope="class")
def prepare_facet_tests():
    if not os.path.isfile(tcf.MODEL_IMAGE):
        wget = f"wget -q http://www.astron.nl/citt/ci_data/wsclean/{tcf.MODEL_IMAGE}"
        check_call(wget.split())

    if not os.path.isfile(tcf.MWA_MOCK_ARCHIVE):
        wget = f"wget -q {tcf.EVERYBEAM_BASE_URL}MWA-single-timeslot.tar.bz2 -O {MWA_MOCK_ARCHIVE}"
        check_call(wget.split())

    os.makedirs(tcf.MWA_MOCK_MS, exist_ok=True)
    check_call(
        f"tar -xf {tcf.MWA_MOCK_ARCHIVE}  -C {tcf.MWA_MOCK_MS} --strip-components=1".split()
    )

    # From python 3.8 onwards, use copytree(..., dirs_exist_ok=True)
    if not os.path.isdir(tcf.MWA_MOCK_FULL):
        shutil.copytree(tcf.MWA_MOCK_MS, tcf.MWA_MOCK_FULL)

    if not os.path.isdir(tcf.MWA_MOCK_FACET):
        shutil.copytree(tcf.MWA_MOCK_MS, tcf.MWA_MOCK_FACET)


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    # Fixture runs only at end of test session, and only cleans-up
    # the working directory if CLEANUP is in your environment variables

    def remove_tarballs():
        for tar in glob.glob("*.tar.*"):
            os.remove(tar)

    def remove_measurement_sets():
        for dir in glob.glob("*.ms"):
            shutil.rmtree(dir)

    def remove_fits_images():
        for f in glob.glob("*.fits"):
            os.remove(f)

    def remove_h5_files():
        for f in glob.glob("*.h5"):
            os.remove(f)

    def collect_cleanup():
        if "CLEANUP" in os.environ:
            os.chdir(tcf.WORKDIR)
            remove_tarballs()
            remove_measurement_sets()
            remove_fits_images()
            remove_h5_files()

    request.addfinalizer(collect_cleanup)


# TODO: remove
# @pytest.fixture(scope="session", autouse=True)
# def cleanup(request):
#     # Fixture is run at end of test session and only does
#     # executes something if CLEANUP to be in your environment
#     def remove_mwa():
#         if "CLEANUP" in os.environ:
#             os.chdir(tcf.WORKDIR)
#             os.remove(MWA_MOCK_ARCHIVE)
#             shutil.rmtree(os.environ["MWA_MOCK_FULL"])
#             shutil.rmtree(os.environ["MWA_MOCK_FACET"])
#             shutil.rmtree(MWA_MOCK_COPY_1, ignore_errors=True)
#             shutil.rmtree(MWA_MOCK_COPY_2, ignore_errors=True)

#     request.addfinalizer(remove_mwa)
