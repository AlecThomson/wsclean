import pytest
import os
import sys
import numpy as np
from utils import validate_call, compute_rms

# Append current directory to system path in order to import testconfig
sys.path.append(".")
# Import configuration variables as test configuration (tcf)
import config_vars as tcf


@pytest.mark.usefixtures("prepare_large_ms")
class TestInlineWeightsAndCorrections:
    def test_inline_apply_weights_and_corrections(self):
        # Test facet-based imaging and applying h5 solutions in an inline manner
        h5download = f"wget -N -q {tcf.WSCLEAN_DATA_URL}/mock_soltab_2pol.h5"
        validate_call(h5download.split())

        name = f"facet-h5-inline-corrections"
        s = f"{tcf.WSCLEAN} -gridder wgridder -name {name} -apply-facet-solutions mock_soltab_2pol.h5 ampl000,phase000 -pol i -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} -niter 1000000 -auto-threshold 5 -mgain 0.8 {tcf.MWA_MS}"
        validate_call(s.split())

        # Check rms of output images is within expected range
        assert os.path.isfile(name + "-dirty.fits")
        rms_dirty = compute_rms(name + "-dirty.fits")
        assert abs(0.358219 - rms_dirty) / rms_dirty < 0.01

        assert os.path.isfile(name + "-model.fits")
        rms_model = compute_rms(name + "-model.fits")
        assert abs(0.10708495 - rms_model) / rms_model < 0.01

        assert os.path.isfile(name + "-residual.fits")
        rms_residual = compute_rms(name + "-residual.fits")
        assert abs(0.26990286 - rms_residual) / rms_residual < 0.01

        assert os.path.isfile(name + "-image.fits")
        rms_image = compute_rms(name + "-image.fits")
        assert abs(0.318041 - rms_image) / rms_image < 0.01
