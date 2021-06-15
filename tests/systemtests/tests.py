import pytest
import os
import warnings
from subprocess import check_call
import shutil

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf

MWA_MS = "MWA-1052736496-averaged.ms"
MWA_COEFFS = ""
CWD = os.getcwd()


@pytest.fixture(autouse=True)
def startup():
    # Fixture is ran for every test function

    # Check if file exist in test_data
    os.makedirs(tcf.WORKDIR, exist_ok=True)
    os.chdir(tcf.WORKDIR)

    if not "MS" in os.environ:
        if not os.path.isfile(f"{MWA_MS}/table.f1"):
            check_call(
                ["wget", "-q", f"www.astron.nl/citt/ci_data/EveryBeam/{MWA_MS}.tgz"]
            )
            check_call(["tar", "-xf", f"{MWA_MS}.tgz"])
            os.remove(f"{MWA_MS}.tgz")
        os.environ["MWA_MS"] = os.path.join(os.getcwd(), MWA_MS)
    # Do the same for (beam) coefficients file

    # change to original working directory
    os.chdir(CWD)


@pytest.fixture(scope="session", autouse=True)
def cleanup(request):
    # Fixture ran only at end of test session

    # Remove measurement set
    def remove_mwa_ms():
        shutil.rmtree(os.environ["MWA_MS"])

    request.addfinalizer(remove_mwa_ms)


def test_dirty_image():
    # Make dirty image
    s = f"wsclean -name test-dirty {tcf.DIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


def test_clean_rectangular_unpadded_image():
    # Clean a rectangular unpadded image
    s = f"wsclean -name clean-rectangular -padding 1 \
          -auto-threshold 5 -mgain 0.8 -niter 100000 {tcf.RECTDIMS} {os.environ['MWA_MS']}"
    check_call(s.split())


# def test_automask_multiscale_clean():
#     # Auto-masked multi-scale clean
#     s = f"wsclean -name multiscale-automasked -auto-threshold 0.5 -auto-mask 3 \
#           -mgain 0.8 -multiscale -niter 100000 {tcf.RECTDIMS} {tcf.MS}"
#     check_call(s.split())

# REQUIRES A LARGER MS
# def test_multiple_intervals():
#     # Multiple intervals
#     s = f"wsclean -name intervals -intervals-out 3 {tcf.RECTDIMS} {tcf.MS}"
#     check_call(s.split())


# def test_multiple_intervals_and_channels():
#     # Multiple intervals + multiple channels with some cleaning
#     s = f"wsclean -name intervals-and-channels -intervals-out 3 \
#         -channels-out 2 -niter 1000 -mgain 0.8 ${tcf.DIMS} {tcf.MS}"
#     check_call(s.split())


# def test_multifrequency_hogbom():
#     # Multi-frequency Högbom clean, no parallel gridding
#     s = f"wsclean -name mfhogbom -channels-out 4 -join-channels -auto-threshold 3 \
#         -mgain 0.8 -niter 1000000 {tcf.RECTDIMS} {tcf.MS}"
#     check_call(s.split())


# def test_multifrequency_hogbom_spectral_fit():
#     # Multi-frequency Högbom clean with spectral fitting
#     s = f"wsclean -name mfhogbom-fitted -channels-out 4 -join-channels -parallel-gridding 4 \
#        -fit-spectral-pol 2 -auto-threshold 3 -mgain 0.8 \
#            -niter 1000000 {tcf.RECTDIMS} {tcf.MS}"
#     check_call(s.split())


# def test_mutifrequency_multiscale_parallel():
#     # Multi-frequency multi-scale clean with spectral fitting, pallel gridding & cleaning
#     s = f"wsclean -name mfms-fitted -channels-out 4 -join-channels -parallel-gridding 4 \
#          -parallel-deconvolution 1000 -fit-spectral-pol 2 -multiscale -auto-threshold 0.5 \
#               -auto-mask 3 -mgain 0.8 -niter 1000000 {tcf.RECTDIMS} {tcf.MS}"
#     check_call(s.split())


# def test_save_components():
#     # Save the list of components
#     s = f"wsclean -name mfms-components -save-source-list -channels-out 4 \
#         -join-channels -parallel-gridding 4 -fit-spectral-pol 2 \
#             -auto-threshold 0.5 -auto-mask 3 -mgain 0.8 -niter 1000000 \
#                 -multiscale -parallel-deconvolution 1000 {tcf.DIMS} {tcf.MS}"
#     check_call(s.split())

# def test_linear_joined_polarizations():
#     # Linear joined polarizations with 4 joined channels
#     s = f"wsclean -name linearpol -niter 1000000 -auto-threshold 3.0 \
#          -pol XX,YY,XY,YX -join-polarizations -join-channels -mgain 0.85 \
#              -channels-out 4 -parallel-gridding 16 {tcf.DIMS} {tcf.MS}"
#     check_call(s.split())

# def test_two_timesteps():
#     # Image two timesteps
#     s = f"wsclean -name two-timesteps -niter 1000000 -auto-threshold 3.0 \
#         -intervals-out 2 -interval 20 22 -mgain 0.85 {tcf.RECTDIMS} {tcf.MS}"
#     check_call(s.split())

# def test_stop_on_negative_components():
#     # Stop on negative components
#     s = f"wsclean -name stop-on-negatives -stop-negative -niter 100000 {tcf.RECTDIMS} {tcf.MS}"
#     check_call(s.split())

# @pytest.mark.parametrize("gridder, name", (["", "shift-ws"], ["-use-wgridder", "shift-wg"]))
# def test_shift_image(gridder, name):
#     # Shift the image with w-stacking and w-gridder gridder
#     s = f"wsclean {gridder} -name {name} -mgain 0.8 -auto-threshold 5 -niter 1000000 -make-psf {tcf.RECTDIMS} -shift 08h09m20s -39d06m54s -no-update-model-required {tcf.MS}"
#     check_call(s.split())

# def test_two_facets():
#     # Apply the facet to the image
#     s = f"wsclean -name two-facets -facet-regions {tcf.FACETFILE_2FACETS} \
#         {tcf.RECTDIMS} {tcf.MS}"
#     check_call(s.split())

# def test_nfacets_pol_xx_yy():
#     # Request two polarizations on approximately 25 facets
#     s = f"wsclean -name nfacets-XX_YY -pol XX,YY \
#         -facet-regions {tcf.FACETFILE_NFACETS} {tcf.RECTDIMS} {tcf.MS}"
#     check_call(s.split())

# def test_facet_beam():
#     if tcf.MWA_PATH:
#         s = f" wsclean -name nfacets-XX_YY-facet-beam -apply-facet-beam -pol XX,YY \
#             -facet-regions {tcf.FACETFILE_NFACETS} {tcf.RECTDIMS} \
#                 -mwa-path {tcf.MWA_PATH} {tcf.MS}"
#         check_call(s.split())
#     else:
#         warnings.warn("MWA_PATH environment variable not set, tests that require beam corrections will be skipped.")


# # MPI
# mpirun wsclean-mp -name mpi ${rectdims} -scale 1amin -channels-out 2 -join-channels -niter 1000000 -mgain 0.8 -auto-threshold 5 -multiscale -no-update-model-required ${ms}
