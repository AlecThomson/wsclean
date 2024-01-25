import pytest
import shutil
import sys
from astropy.io import fits
from astropy.wcs import WCS
import casacore.tables
import h5py
import numpy as np
from utils import (
    assert_taql,
    basic_image_check,
    check_and_remove_files,
    compare_rms_fits,
    validate_call,
)

# Append current directory to system path in order to import testconfig
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


def gridders():
    return {
        "wstacking": "-gridder wstacking",
        "wgridder": "-gridder wgridder",
        "idg": "-gridder idg",
    }


def predict_full_image(ms, gridder):
    """Predict full image"""
    s = f"{tcf.WSCLEAN} -predict -gridder {gridder} -name point-source {ms}"
    validate_call(s.split())


def predict_facet_image(ms, gridder, apply_facet_beam):
    name = "point-source"
    facet_beam = "-apply-facet-beam -mwa-path ." if apply_facet_beam else ""
    if apply_facet_beam:
        shutil.copyfile(f"{name}-model.fits", f"{name}-model-pb.fits")

    # Predict facet based image
    s = f"{tcf.WSCLEAN} -predict -gridder {gridder} {facet_beam} -facet-regions {tcf.FACETFILE_4FACETS} -name {name} {ms}"
    validate_call(s.split())


def deconvolve_facets(ms, gridder, reorder, mpi, apply_beam=False):
    nthreads = 4
    mpi_cmd = f"mpirun -tag-output -np {nthreads} {tcf.WSCLEAN_MP}"
    thread_cmd = f"{tcf.WSCLEAN} -parallel-gridding {nthreads}"
    reorder_ms = "-reorder" if reorder else "-no-reorder"
    s = [
        mpi_cmd if mpi else thread_cmd,
        f"-gridder {gridder}",
        tcf.DIMS_SMALL,
        reorder_ms,
        "-niter 1000000 -auto-threshold 5 -mgain 0.8",
        f"-facet-regions {tcf.FACETFILE_4FACETS}",
        f"-name facet-imaging{reorder_ms}",
        "-mwa-path . -apply-facet-beam" if apply_beam else "",
        "-v",
        ms,
    ]
    validate_call(" ".join(s).split())


def create_pointsource_grid_skymodel(
    skymodel_filename, grid_size, nr_pixels, wcs
):
    """
    Writes a skymodel file for a square grid of point sources in a square image.

    Parameters
    ----------
    skymodel_filename: str
    grid_size: int
        Number of point sources (in one direction)
    nr_pixels: int
        Number of pixels in the image (in one direction)
    wcs: astropy.wcs.WCS
        World coordinate system of the image

    Returns
    -------
    list of tuples
        A list of source/pixel positions of length grid_size*grid_size
    """
    source_positions = []
    source_pixel_index_range = (np.arange(grid_size)) * (
        nr_pixels // grid_size
    ) + (nr_pixels // grid_size // 2)
    with open(skymodel_filename, "w") as sky_model_file:
        print(
            "Format = Name, Patch, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='150000000', MajorAxis, MinorAxis, Orientation",
            file=sky_model_file,
        )
        for i, idx0 in enumerate(source_pixel_index_range):
            for j, idx1 in enumerate(source_pixel_index_range):
                sky = wcs.pixel_to_world(idx0, idx1, 0, 0)
                print(
                    f",direction_{i}{j},,{sky[0].ra.rad},{sky[0].dec.rad},,,,,,,",
                    file=sky_model_file,
                )
                print(
                    f"source-{i}-{j},direction_{i}{j},POINT,{sky[0].ra.rad},{sky[0].dec.rad},1.0,[],false,150000000,,,",
                    file=sky_model_file,
                )
                source_positions.append((idx0, idx1))

    return source_positions


@pytest.mark.usefixtures("prepare_mock_ms", "prepare_model_image")
class TestFacets:
    def test_makepsfonly(self):
        """
        Test that wsclean with the -make-psf-only flag exits gracefully and
        that the psf passes basic checks.
        """
        s = [
            tcf.WSCLEAN,
            "-make-psf-only -name facet-psf-only",
            tcf.DIMS_SMALL,
            f"-facet-regions {tcf.FACETFILE_4FACETS} {tcf.MWA_MOCK_MS}",
        ]
        validate_call(" ".join(s).split())

        basic_image_check("facet-psf-only-psf.fits")

    # Test assumes that IDG and EveryBeam are installed
    @pytest.mark.parametrize("gridder", gridders().items())
    def test_stitching(self, gridder):
        """Test stitching of the facets"""
        prefix = f"facet-stitch-{gridder[0]}"
        s = [
            tcf.WSCLEAN,
            "-quiet",
            gridder[1],
            tcf.DIMS_SMALL,
            "" if (gridder[0] == "idg") else "-pol XX,YY",
            f"-facet-regions {tcf.FACETFILE_2FACETS}",
            f"-name {prefix}",
            tcf.MWA_MOCK_MS,
        ]
        validate_call(" ".join(s).split())
        fpaths = (
            [prefix + "-dirty.fits", prefix + "-image.fits"]
            if (gridder[0] == "idg")
            else [
                prefix + "-XX-dirty.fits",
                prefix + "-YY-dirty.fits",
                prefix + "-XX-image.fits",
                prefix + "-YY-image.fits",
            ]
        )
        check_and_remove_files(fpaths, remove=True)

    # FIXME: we should test wstacking and wgridder here too
    # but they fail on the taql assertion
    @pytest.mark.parametrize("gridder", ["wgridder"])
    @pytest.mark.parametrize("apply_facet_beam", [False, True])
    def test_predict(self, gridder, apply_facet_beam):
        """
        Test predict only run

        Parameters
        ----------
        gridder : str
            wsclean compatible description of gridder to be used.
        """

        predict_facet_image(tcf.MWA_MOCK_FACET, gridder, apply_facet_beam)

        # A numerical check can only be performed in case no DD effects were applied.
        if not apply_facet_beam:
            predict_full_image(tcf.MWA_MOCK_FULL, gridder)
            taql_command = f"select from {tcf.MWA_MOCK_FULL} t1, {tcf.MWA_MOCK_FACET} t2 where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,5e-3))"
            assert_taql(taql_command)

    @pytest.mark.parametrize("gridder", ["wgridder"])
    @pytest.mark.parametrize("reorder", [False, True])
    @pytest.mark.parametrize("mpi", [False, True])
    def test_facetdeconvolution(self, gridder, reorder, mpi):
        """
        Test facet-based deconvolution

        Parameters
        ----------
        gridder : str
            wsclean compatible description of gridder to be used.
        reorder : bool
            Reorder MS?
        mpi : bool
            True: Use MPI for parallel gridding.
            False: Use multi-threading for parallel gridding.
        """
        # Parametrization causes some overhead in that predict of full image is run for
        # every parametrization
        predict_full_image(tcf.MWA_MOCK_FULL, gridder)

        # Make sure old versions of the facet mock ms are removed
        shutil.rmtree(tcf.MWA_MOCK_FACET)

        # Copy output to new MS, swap DATA column, and remove MODEL_DATA
        validate_call(
            f"cp -r {tcf.MWA_MOCK_FULL} {tcf.MWA_MOCK_FACET}".split()
        )
        assert shutil.which("taql") is not None, "taql executable not found!"

        validate_call(
            [
                "taql",
                "-noph",
                f"UPDATE {tcf.MWA_MOCK_FACET} SET DATA=MODEL_DATA",
            ]
        )
        validate_call(
            [
                "taql",
                "-noph",
                f"ALTER TABLE {tcf.MWA_MOCK_FACET} DROP COLUMN MODEL_DATA",
            ]
        )
        taql_command = f"select from {tcf.MWA_MOCK_FULL} t1, {tcf.MWA_MOCK_FACET} t2 where not all(near(t1.MODEL_DATA,t2.DATA, 4e-3))"
        assert_taql(taql_command)

        deconvolve_facets(tcf.MWA_MOCK_FACET, gridder, reorder, mpi)

        taql_command = f"select from {tcf.MWA_MOCK_FACET} where not all(near(DATA,MODEL_DATA, 4e-3))"
        assert_taql(taql_command)

    def test_read_only_ms(self):
        chmod = f"chmod a-w -R {tcf.MWA_MOCK_FULL}"
        validate_call(chmod.split())
        try:
            # When "-no-update-model-required" is specified, processing a read-only measurement set should be possible.
            s = f"{tcf.WSCLEAN} -interval 10 20 -no-update-model-required -name facet-readonly-ms -auto-threshold 0.5 -auto-mask 3 \
                -mgain 0.95 -nmiter 2 -multiscale -niter 100000 \
                -facet-regions {tcf.FACETFILE_4FACETS} \
                {tcf.DIMS_SMALL} {tcf.MWA_MOCK_FULL}"
            validate_call(s.split())
        finally:
            chmod = f"chmod u+w -R {tcf.MWA_MOCK_FULL}"
            validate_call(chmod.split())

    @pytest.mark.parametrize("mpi", [False, True])
    def test_facetbeamimages(self, mpi):
        """
        Basic checks of the generated images when using facet beams. For each image,
        test that the pixel values are valid (not NaN/Inf) and check the percentage
        of zero pixels.
        """

        deconvolve_facets(tcf.MWA_MOCK_FACET, "wgridder", True, mpi, True)

        basic_image_check("facet-imaging-reorder-psf.fits")
        basic_image_check("facet-imaging-reorder-dirty.fits")

    def test_multi_channel(self):
        # Test for issue 122. Only test if no crash occurs.
        h5download = f"wget -N -q {tcf.WSCLEAN_DATA_URL}/mock_soltab_2pol.h5"
        validate_call(h5download.split())

        s = f"{tcf.WSCLEAN} -parallel-gridding 3 -channels-out 2 -gridder wgridder -name multi-channel-faceting -apply-facet-solutions mock_soltab_2pol.h5 ampl000,phase000 -pol xx,yy -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} -join-polarizations -interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 {tcf.MWA_MOCK_MS}"
        validate_call(s.split())

    def test_diagonal_solutions(self):
        h5download = f"wget -N -q {tcf.WSCLEAN_DATA_URL}/mock_soltab_2pol.h5"
        validate_call(h5download.split())

        s = f"{tcf.WSCLEAN} -parallel-gridding 3 -channels-out 2 -gridder wgridder -name faceted-diagonal-solutions -apply-facet-solutions mock_soltab_2pol.h5 ampl000,phase000 -diagonal-solutions -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} -interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 {tcf.MWA_MOCK_MS}"
        validate_call(s.split())

    def test_diagonal_solutions_with_beam(self):
        h5download = f"wget -N -q {tcf.WSCLEAN_DATA_URL}/mock_soltab_2pol.h5"
        validate_call(h5download.split())

        s = f"{tcf.WSCLEAN} -parallel-gridding 3 -channels-out 2 -gridder wgridder -name faceted-diagonal-solutions -apply-facet-solutions mock_soltab_2pol.h5 ampl000,phase000 -diagonal-solutions -mwa-path . -apply-facet-beam -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} -interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 {tcf.MWA_MOCK_MS}"
        validate_call(s.split())

    def test_parallel_gridding(self):
        # Compare serial, threaded and mpi run for facet based imaging
        # with h5 corrections. Number of used threads/processes is
        # deliberately chosen smaller than the number of facets.
        h5download = f"wget -N -q {tcf.WSCLEAN_DATA_URL}/mock_soltab_2pol.h5"
        validate_call(h5download.split())

        names = ["facets-h5-serial", "facets-h5-threaded", "facets-h5-mpi"]
        wsclean_commands = [
            tcf.WSCLEAN,
            f"{tcf.WSCLEAN} -parallel-gridding 3",
            f"mpirun -np 3 {tcf.WSCLEAN_MP}",
        ]
        for name, command in zip(names, wsclean_commands):
            s = f"{command} -name {name} -apply-facet-solutions mock_soltab_2pol.h5 ampl000,phase000 -pol xx,yy -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} -join-polarizations -interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 -gridder wstacking {tcf.MWA_MOCK_MS}"
            validate_call(s.split())

        # Typical rms difference is about 1.0e-7
        threshold = 3.0e-7
        compare_rms_fits(
            f"{names[0]}-YY-image.fits", f"{names[1]}-YY-image.fits", threshold
        )
        compare_rms_fits(
            f"{names[0]}-YY-image.fits", f"{names[2]}-YY-image.fits", threshold
        )

    @pytest.mark.parametrize("beam", [False, True])
    @pytest.mark.parametrize(
        "h5file",
        [
            None,
            ["mock_soltab_2pol.h5"],
            ["mock_soltab_2pol.h5", "mock_soltab_2pol.h5"],
        ],
    )
    def test_multi_ms(self, beam, h5file):
        """
        Check that identical images are obtained in case multiple (identical) MSets and H5Parm
        files are provided compared to imaging one MSet
        """

        h5download = f"wget -N -q {tcf.WSCLEAN_DATA_URL}/mock_soltab_2pol.h5"
        validate_call(h5download.split())

        # Make a new copy of tcf.MWA_MOCK_MS into two MSets
        validate_call(f"cp -r {tcf.MWA_MOCK_MS} {tcf.MWA_MOCK_COPY_1}".split())
        validate_call(f"cp -r {tcf.MWA_MOCK_MS} {tcf.MWA_MOCK_COPY_2}".split())

        names = ["facets-single-ms", "facets-multiple-ms"]
        commands = [
            f"{tcf.MWA_MOCK_MS}",
            f"{tcf.MWA_MOCK_COPY_1} {tcf.MWA_MOCK_COPY_2}",
        ]

        if beam:
            commands = [
                "-mwa-path . -apply-facet-beam " + command
                for command in commands
            ]

        if h5file is not None:
            commands[0] = (
                f"-apply-facet-solutions {h5file[0]} ampl000,phase000 "
                + commands[0]
            )
            commands[1] = (
                f"-apply-facet-solutions {','.join(h5file)} ampl000,phase000 "
                + commands[1]
            )

        # Note: -j 1 enabled to ensure deterministic iteration over visibilities
        for name, command in zip(names, commands):
            s = f"{tcf.WSCLEAN} -j 1 -nmiter 2 -gridder wgridder -name {name} -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} -interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 {command}"
            validate_call(s.split())

        # Compare images.
        threshold = 1.0e-6
        compare_rms_fits(
            f"{names[0]}-image.fits", f"{names[1]}-image.fits", threshold
        )

        # Model data columns should be equal
        taql_commands = [
            f"select from {tcf.MWA_MOCK_MS} t1, {tcf.MWA_MOCK_COPY_1} t2 where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,1e-6))"
        ]
        taql_commands.append(
            f"select from {tcf.MWA_MOCK_COPY_1} t1, {tcf.MWA_MOCK_COPY_2} t2 where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,1e-6))"
        )
        # assert_taql(taql_command for taql_command in taql_commands)
        for taql_command in taql_commands:
            assert_taql(taql_command)

    def test_diagonal_solutions(self):
        # Initialize random rumber generator
        rng = np.random.default_rng(1)

        # Strip unused stations from mock measurement set
        s = f"DP3 msin={tcf.MWA_MOCK_MS}  msout=diagonal_solutions.ms msout.overwrite=True steps=[filter] filter.remove=True"
        validate_call(s.split())

        # Fill WEIGHT_SPECTRUM with random values
        with casacore.tables.table(
            "diagonal_solutions.ms", readonly=False
        ) as t:
            weight_spectrum_shape = np.concatenate(
                (
                    np.array([t.nrows()]),
                    t.getcoldesc("WEIGHT_SPECTRUM")["shape"],
                )
            )
            weights = rng.uniform(0, 1, weight_spectrum_shape) + np.array(
                [1, 2, 3, 4], ndmin=3
            )
            t.putcol("WEIGHT_SPECTRUM", weights)

        # Create a template image
        s = (
            f"{tcf.WSCLEAN} -gridder wgridder -name template-diagonal-solutions "
            f"-size 256 256 -scale 4amin -interval 0 1 diagonal_solutions.ms"
        )
        validate_call(s.split())

        # Use template image to create a sky model consisting of a grid of point sources
        with fits.open("template-diagonal-solutions-image.fits") as f:
            wcs = WCS(f[0].header)
            nr_pixels = f[0].shape[-1]
        pointsource_grid_size = 2
        source_positions = create_pointsource_grid_skymodel(
            "diagonal-solutions-skymodel.txt",
            pointsource_grid_size,
            nr_pixels,
            wcs,
        )

        # Predict (without solutions)
        s = f"DP3 msin=diagonal_solutions.ms msout= steps=[predict] predict.sourcedb=diagonal-solutions-skymodel.txt"
        validate_call(s.split())

        # Image (without solutions)
        s = (
            f"{tcf.WSCLEAN} -name diagonal-solutions-reference -size 256 256 "
            "-no-reorder "
            "-scale 4amin "
            "diagonal_solutions.ms"
        )
        validate_call(s.split())

        # Create template solutions .h5 file
        s = "DP3 msin=diagonal_solutions.ms msout= steps=[ddecal] ddecal.sourcedb=diagonal-solutions-skymodel.txt ddecal.h5parm=diagonal-solutions.h5 ddecal.mode=complexgain"
        validate_call(s.split())

        # Fill the template solutions file with random data
        with h5py.File("diagonal-solutions.h5", mode="r+") as f:
            f["sol000"]["phase000"]["val"][:] = rng.uniform(
                -np.pi, np.pi, f["sol000"]["phase000"]["val"].shape
            )
            f["sol000"]["phase000"]["weight"][:] = 1.0
            f["sol000"]["amplitude000"]["val"][:] = rng.uniform(
                0.5, 3, f["sol000"]["amplitude000"]["val"].shape
            )
            f["sol000"]["amplitude000"]["weight"][:] = 1.0

        # Predict with (random) solutions
        s = (
            "DP3 msin=diagonal_solutions.ms msout= steps=[h5parmpredict] "
            "h5parmpredict.sourcedb=diagonal-solutions-skymodel.txt "
            "h5parmpredict.applycal.parmdb=diagonal-solutions.h5 "
            "h5parmpredict.applycal.steps=[ampl,phase] "
            "h5parmpredict.applycal.ampl.correction=amplitude000 "
            "h5parmpredict.applycal.phase.correction=phase000 "
            "h5parmpredict.applycal.correction=amplitude000"
        )
        validate_call(s.split())

        # Image data predicted with solutions applied,
        # without applying corrections for the solutions while imaging
        s = (
            f"{tcf.WSCLEAN} -name diagonal-solutions-no-correction -size 256 256 "
            "-no-reorder "
            "-scale 4amin "
            "diagonal_solutions.ms"
        )
        validate_call(s.split())

        # Image data predicted with solutions applied,
        # while applying corrections
        s = (
            f"{tcf.WSCLEAN} -name diagonal-solutions "
            "-parallel-gridding 3 "
            "-size 256 256 "
            "-no-reorder "
            "-scale 4amin "
            "-mgain 0.8 "
            "-threshold 10mJy "
            "-niter 10000 "
            f"-facet-regions {tcf.FACETFILE_4FACETS} "
            "-apply-facet-solutions diagonal-solutions.h5 "
            "amplitude000,phase000 -diagonal-solutions "
            "diagonal_solutions.ms"
        )
        validate_call(s.split())

        # Compare reference, uncorrection and corrected fluxes
        reference_image_data = fits.getdata(
            "diagonal-solutions-reference-image.fits"
        )[0, 0]
        no_correction_image_data = fits.getdata(
            "diagonal-solutions-no-correction-image.fits"
        )[0, 0]
        image_data = fits.getdata("diagonal-solutions-image-pb.fits")[0, 0]
        # loop over input sources
        for idx0, idx1 in source_positions:
            # Assert that without corrections less than 5 percent flux is recovered
            assert np.abs(no_correction_image_data[idx0, idx1]) < 5e-2
            # Assert that with corrections the recovered flux is within 2 percent of the reference
            assert np.isclose(
                reference_image_data[idx0, idx1],
                image_data[idx0, idx1],
                rtol=2e-2,
            )

    def test_dd_psfs_with_faceting(self):
        s = [
            tcf.WSCLEAN,
            "-name dd-psfs-with-faceting -dd-psf-grid 3 3 -parallel-gridding 5 -parallel-deconvolution 100 -channels-out 2 -join-channels -niter 100 -mgain 0.8 -gridder wgridder -apply-facet-beam -mwa-path .",
            tcf.DIMS_SMALL,
            f"-facet-regions {tcf.FACETFILE_4FACETS} {tcf.MWA_MOCK_MS}",
        ]
        validate_call(" ".join(s).split())
        import os.path

        basic_image_check("dd-psfs-with-faceting-MFS-image.fits")
        for i in range(9):
            assert os.path.isfile(
                f"dd-psfs-with-faceting-d000{i}-0000-psf.fits"
            )
            assert os.path.isfile(
                f"dd-psfs-with-faceting-d000{i}-0001-psf.fits"
            )
            assert os.path.isfile(
                f"dd-psfs-with-faceting-d000{i}-MFS-psf.fits"
            )
        assert not os.path.isfile(f"dd-psfs-with-faceting-0000-psf.fits")
        assert not os.path.isfile(f"dd-psfs-with-faceting-0001-psf.fits")
        assert not os.path.isfile(f"dd-psfs-with-faceting-MFS-psf.fits")
