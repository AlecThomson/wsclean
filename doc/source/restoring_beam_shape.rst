Restoring beam shape
====================

This chapter describes WSClean's algorithm for determining the synthesized beam major and minor axes and the position angle (PA) that is used for restoring, and presents some solutions if the restoring beam does not have a desirable size or shape.

Solving highly elliptical beams
-------------------------------

When the fit is highly elliptical and when this is undesirable, the option ``-circular-beam`` can be added. This constrains the fit to be circular. When forcing a circular beam, WSClean makes still sure that the resulting flux scale of restored sources is correct. This may improve the esthetics of produced images, but it might of course give a false sense of having a circular beam.

If neither the default nor a circular fit is adequate, another option is to specify the beam size (=just a radius) or even shape (major and minor axis, and position angle) manually. WSClean provides options ``-beam-size`` and ``-beam-shape`` to set the size or shape manually. Two examples:

.. code-block:: bash

    wsclean -beam-size 1.5amin [options...] observation.ms

This will restore the model on top of the residual using a circular Gaussian of size 1.5'.

.. code-block:: bash

    wsclean -beam-shape 30asec 20asec 45deg [options...] observation.ms

This will restore using an elliptical Gaussian with major size of 30", minor size of 20", rotated by 45 degree.

Fitting algorithm
-----------------

The input to the fitting algorithm is the PSF image and an estimate of the PSF FWHM that is determined from the longest baseline size and frequency that is gridded. Given this first estimate of the beamsize, a box width and height is calculated to be :math:`10 \times` the estimated FWHM of the PSF around the central pixel. A 3-parameter least-squares fit is then performed to find the beam's major, minor and PA values inside that box, to minimize the sum-of-squares difference between the fitted model and the real PSF over all pixels in that box of :math:`10 \times` the size of the PSF. So this includes values out to a few lobes.

Using the full image to do the fit is too expensive — the factor of 10 is to make the fitting go fast. The fit is generally equal to a fit using a larger area.

There are some added heuristical things when fits are very different from the initial estimated PSF: WSClean will iteratively recalculate the fitting box and refit.

History
-------

* In :doc:`WSClean 2.6 <changelogs/v2.6>`, the option ``-beam-fitting-size`` was added. It sets the box size to be used during the fit. If the beam that is fitted is different from what is expected, changing this value might help. It seems that in particular lowering the value to 1-3 can solve certain (rare) issues.
* In :doc:`WSClean 2.5 <changelogs/v2.5>`, an improved fit for small and forced circular beams was added.
