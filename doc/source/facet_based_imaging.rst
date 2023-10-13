Facet-based imaging
===================

WSClean supports facet-based imaging with correction of direction-dependent effects.
A facet is a (polygonal-shaped) subsection of the full image. WSClean supports both convex and concave polygons.

Facet-based imaging supports the correction of each facet with a previously found (direction-dependent) gain solution. 
In some cases, facet-based imaging can also speed up imaging or reduce the memory footprint.

Command-line options
--------------------

To enable facet-based imaging in WSClean, a file containing the facet definitions should be provided via the :code:`-facet-regions` option on the command line:

.. code-block:: text

    -facet-regions <MY_REGIONS_FILE>

in which :code:`[MY_REGIONS_FILE]` is a file containing the facet definitions in the DS9 region file format.
For a detailed explanation on the expected file format, see the explanation on :doc:`ds9_facet_file`.

Enabling the facet beam correction can be done with the option

.. code-block:: text

    -apply-facet-beam

The facet beam update interval (in seconds) can be defined by specifying:

.. code-block:: text

    -facet-beam-update <seconds>

The default value for the update interval is 120s.

Direction-dependent corrections per facet can be read and applied from an ``h5parm`` file, which is a HDF5 file with gain solutions in a particular format. This format is for example supported by tools like `DP3 <https://dp3.readthedocs.io/>`_ and `losoto <https://github.com/revoltek/losoto>`_.

The ``h5parm`` file is specified via the command-line option:

.. code-block:: text

    -apply-facet-solutions <path-to-h5parm> <name1[,name2]>

where :code:`<path-to-h5parm>` is the path to the ``h5parm`` solution file and :code:`<name1[,name2]>`
is a comma-separated list of strings specifying which "soltabs" from the provided ``h5parm`` solution file are used.
Acceptable names are :code:`ampl000` and/or :code:`phase000`.

WSClean supports application of ``h5parm`` with various polarizations:

- The simplest way is to apply "scalar" solutions. For example, if the solution file specifies ``X`` or ``R`` solutions, a scalar correction can be applied to create corrected ``XX`` or ``RR`` images, respectively. If the solution file specifies multiple polarizations, it is possible to image and correct multiple polarization at once (e.g. ``-pol xx,yy``). This way, diagonal solutions can be applied, but it requires imaging the ``xx`` and ``yy`` images separately, which is not always desired.
- Stokes I images with diagonal solutions can be made without imaging ``xx`` and ``yy`` separately by using the option ``-diagonal-solutions``. This makes WSClean read the ``xx`` and ``yy`` visibilities (or ``rr`` and ``ll``), apply the corresponding solutions to them and average them together before gridding.
- Full-polarization IQUV imaging with diagonal or full-Jones solutions is still a work in progress.

In case multiple measurement sets are specified, it is possible to either specify one ``h5parm`` solution file, or a separate ``h5parm`` solution
file per seasurement set.
The correction that should be applied (:code:`ampl000`, :code:`phase000`, or both) is required to be identical for all ``h5parm`` solution files.
As an illustration, assume that :code:`N` measurement sets are passed to WSClean, with corresponding solution files :code:`h5parm1.h5, h5parm2.h5, ..., h5parmN.h5` containing a
scalar amplitude correction.
The syntax for applying the facet solution files on its corresponding measurement set thus becomes:

.. code-block:: text

    -apply-facet-solutions h5parm1.h5,h5parm2.h5,...,h5parmN.h5 ampl000

.. note::
    To find the matching direction in the solution file for the specified facets,
    the (RA, Dec) pointing of each facet is matched against the direction with
    the smallest (Euclidean) distance in the solution file.
    For further information on the (RA, Dec) pointing of a facet, see :doc:`ds9_facet_file`.


Examples
--------
This is an example facet-based imaging command that applies both a facet-based beam correction and a scalar gain correction from an ``h5parm`` file:

.. code-block:: bash

    wsclean \
    -apply-facet-solutions mock_soltab_2pol.h5 ampl000,phase000 \
    -facet-regions ds9.reg \
    -apply-facet-beam \
    -facet-beam-update 120 \
    -niter 1000000 -auto-threshold 5 -mgain 0.8 \
    -size 1024 1024 -scale 1amin \
    ${ms}

In case the solution files contains separate ``x`` and ``y`` solutions, option ``-diagional-solutions`` should be added.
    
Availability
------------
Initial support for facetting is made available in WSClean :doc:`version 3.0 <changelogs/v3.0>`. In subsequent versions,
several bugs were fixed and support for different solution types was added. WSClean :doc:`version 3.4 <changelogs/v3.4>`
has support for scalar and diagonal solutions, and is considered stable.

Facet-based imaging in conjunction with the Image Domain Gridder (IDG) is only possible without applying DDEs.
