import sys
from pathlib import Path

import numpy as np
import pytest
from utils import compute_rms, validate_call

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

        name = "facet-h5-inline-corrections"
        command = f"{tcf.WSCLEAN} -gridder wgridder -name {name} -apply-facet-solutions mock_soltab_2pol.h5 ampl000,phase000 -pol i -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} -niter 1000000 -auto-threshold 5 -mgain 0.8 {tcf.MWA_MS}"
        validate_call(command.split())

        # Check rms of output images is within 1% of the expected rms calculated from previous runs
        expected_rms = {
            "dirty": 0.358219,
            "model": 0.10708495,
            "residual": 0.26990286,
            "image": 0.318041,
        }
        for image_type, expected_image_rms in expected_rms.items():
            file_path = Path(f"{name}-{image_type}.fits")
            assert file_path.is_file()
            rms = compute_rms(file_path)
            assert abs(expected_image_rms - rms) / rms < 0.01
