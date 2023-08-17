Basic cleaning
==============

Image dimensions
----------------

The image size is set in pixels with the ``-size`` parameter, which takes a width and a height. An image is not required to be square shaped. It is however required to use even numbers for width and height. The ``-scale`` parameter takes an angle and sets the angular size of a single pixel. The pixel scale currently has to be square shaped. Together, the ``-size`` and ``-scale`` parameters set the angular size of the image. 

Stopping criteria
-----------------

To enable cleaning, the ``-niter`` parameter should be set to a non-zero value. This parameter sets the maximum number of minor iterations that are allowed to be used. Once the specified number of iterations are reached, cleaning will stop.

It is good practice to let cleaning reach a threshold, and only use ``-niter`` to make sure wsclean will not run for an excessively long time. One should also not clean deeper than the noise, unless a mask is used. The easiest way of setting a stopping criterium based on the noise, is by using an automatic threshold. This is discussed next.

Automatic threshold
-------------------

The recommended way to set a stopping criterion, is by setting a threshold relative to the residual noise level. The option for this is ``-auto-threshold``. With this option, WSClean will calculate the standard deviation of the residual image before the start of every major deconvolution iteration, and clean until the given auto-threshold factor times the standard deviation of the image. The standard deviation is estimated using the medium absolute deviation, which is a robust estimator that is not very sensitive to source structure still present in the image. When performing :doc:`wideband <wideband_deconvolution>` and/or :doc:`polarized deconvolution <polarimetric_deconvolution>`, the standard deviation is measured from the integrated image. An example:

.. code-block:: bash

    # Clean until reaching a 3 sigma noise level
    wsclean -auto-threshold 3 -size 2048 2048 -scale 1amin \
      -mgain 0.8 -niter 50000 observation.ms

Note that the ``-mgain`` parameter is used, in order to enable the Cotton-Schwab style major iterations. While this is not necessary for the automatic threshold to work, if the image has a very high dynamic range, the initially computed standard deviation might not be a good estimate. By using Cotton-Schwab, the standard deviation is recalculated at the beginning of every major iteration, and this will be more accurate. The ``-mgain`` parameter is discussed in more detail in the next section.

In simple runs, a typical stopping criterion is 3 x standard deviation of the noise in the image, and ``-auto-threshold 3`` is therefore a common setting. A more advanced and deeper cleaning can be performed using :doc:`auto-masking <masking>`.

Manual threshold
----------------

The auto-threshold sets a value that is relative to the noise. This is generally preferred over specifying an absolute threshold flux-density, because it avoids having to know the level of the noise before running the clean. An automated procedure for chosing the threshold also avoids bias and improves reproducability. 

Nevertheless, in some situations it might be desirable to use a manual absolute threshold. For this, the ``-abs-threshold <value>`` parameter can be used. If no unit is specified, the value is specified in units of Janskys. Specifying ``-abs-threshold 0.7`` results thus in a 700 mJy absolute threshold. The unit can also me explicitly specified, e.g. ``-abs-threshold 700mJy``. 

It is possible to specify both an automatic threshold and a manual threshold. In this case, whenever one of the thresholds is reached, cleaning stops.

Using Cotton-Schwab: the ``-mgain`` parameter
---------------------------------------------

The ``-mgain`` parameter sets the major iteration gain: during every major iteration, the peak is reduced by the given factor. With an ``mgain`` of 0.8, the peak is reduced by 80%. This is quite a common and safe option for cleaning. With a reasonable good PSF, using 0.9 is possible without loss of accuracy, and a bit faster. With a very bad PSF, it might be necessary to lower the ``mgain`` parameter.

When ``-mgain`` is not given, or when it is set to 1, WSClean will never go back to the visibilities. It will therefore perform a simple image-based Högbom clean. While this is fast, it limits the accuracy of the image and the dynamic range that can be reached. WSClean will also not write to the MODEL column, as would be required for self-calibration (see the chapter on :doc:`self-calibration <selfcal>`). To use Högbom clean and still fill the model column, an mgain of e.g. 0.9999 can be used. Alternatively, the final model :doc:`can be predicted <prediction>` in a second WSClean run.

Note that the ``-mgain`` parameter is not the same as the ``-gain``  parameter. The latter sets the minor loop cleaning gain, and is 0.1 by default. It is almost never required to change the ``gain`` parameter.

**Next chapter:** :doc:`prediction`
