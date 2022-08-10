#!/bin/sh
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# Author: Chiara Salvoni

# Script for downloading a 9h long SKA-mid simulated measurement set.
# This set contains a point source off the phase center, and it is heavily averaged.

set -e

# Download into a 'test_data' directory inside the current directory.
mkdir -p test_data
cd test_data/

SKA_LOW_MS_ARCHIVE=ska-low-sim-averaged.tar.gz
SKA_LOW_MS=ska-low-sim-averaged.ms

if [ ! -f ${SKA_LOW_MS}/table.f0 ] ; then

    if [ ! -f "$SCREEN_FITTING_MS_ARCHIVE" ]; then
        wget -q https://support.astron.nl/software/ci_data/EveryBeam/ska-low-sim-averaged.tar.gz -O $SKA_LOW_MS_ARCHIVE
    fi

    rm -rf $SKA_LOW_MS
    tar xf $SKA_LOW_MS_ARCHIVE
    rm $SKA_LOW_MS_ARCHIVE
fi

