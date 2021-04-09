#!/bin/bash

# Tests the by taking the following sequence of steps:
# - Predict the visibilities from provided model image (without facetting) and write predicted visibilities to MODEL_DATA column of MS_0
# - Copy MS_0 to MS_1 that will be used for facet based imaging
# - Run an imaging cycle (including a number of major iterations) using facets
# - Check whether MODEL_DATA column of MS_0 and MS_1 are close

set -e

# Common settings
facetfile="ds9_facet.reg"

predict_settings="-predict -use-wgridder -name point-source"
rectdims="-size 256 256 -scale 4amin"
major_cycle="-niter 1000000 -auto-threshold 5 -mgain 0.8"

# Make sure that wsclean executable from build directory is used
export PATH=${PWD}:$PATH
cd ${PWD}/test_data

MWA_MOCK_ARCHIVE=MWA_ARCHIVE.tar.bz2
MWA_MOCK_MS=MWA_MOCK.ms
# MS used for full image and facet based predict respectively
MWA_MOCK_FULL=MWA_MOCK_FULL.ms
MWA_MOCK_FACET=MWA_MOCK_FACET.ms

# Download a mock measurement set from the citt server
if [ ! -f ${MWA_MOCK_MS}/table.f1 ]; then

    if [ ! -f "$MWA_MOCK_ARCHIVE" ]; then
	wget -q www.astron.nl/citt/EveryBeam/MWA-single-timeslot.tar.bz2 -O $MWA_MOCK_ARCHIVE
    fi

    mkdir -p $MWA_MOCK_MS

    tar -xf $MWA_MOCK_ARCHIVE  -C $MWA_MOCK_MS --strip-components=1
    rm $MWA_MOCK_ARCHIVE
fi

cp -r $MWA_MOCK_MS $MWA_MOCK_FULL

echo "===== Predicting facet image ====="
wsclean ${predict_settings} ${MWA_MOCK_FULL}

echo "===== Copy to new MeasurementSet ====="
DPPP msin=${MWA_MOCK_FULL} msin.datacolumn="MODEL_DATA" msout=${MWA_MOCK_FACET} msout.datacolumn="DATA" msout.overwrite=true steps=[]

echo "===== Copy to new MeasurementSet ====="
# TODO: remove following line
# wsclean -use-wgridder ${major_cycle} -name facet-imaging ${rectdims} ${MWA_MOCK_FACET}
wsclean -use-wgridder ${major_cycle} -facet-regions ${facetfile} -name facet-imaging ${rectdims} ${MWA_MOCK_FACET}

# taql command to compare
taql "select from MWA_MOCK_FACET.ms t1, MWA_MOCK_FACET.ms t2 where not all(near(t1.DATA,t2.MODEL_DATA, 4e-3))" > taql.out

# Create expected taql output.
echo "    select result of 0 rows" > taql.ref
diff taql.out taql.ref  ||  exit 1

rm -rf ${MWA_MOCK_FULL} ${MWA_MOCK_FACET}