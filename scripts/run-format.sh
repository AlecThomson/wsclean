#!/bin/bash

# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

#Script configuration for this repo. Adjust it when copying to a different repo.

# Disable globbing. This is needed when defining patterns that have wildcards.
set -e -f

#The directory that contains the source files.
SOURCE_DIR=$(dirname "$0")/..

#Directories that must be excluded from formatting. These paths are
#relative to SOURCE_DIR.
EXCLUDE_DIRS=(external build CMake)

#The patterns of the C++ source files, which clang-format should format.
CXX_SOURCES=(*.cpp *.h)

#The patterns of the CMake source files, which cmake-format should format.
CMAKE_SOURCES=(CMakeLists.txt *.cmake)

#End script configuration.

#The common formatting script has further documentation.
source $(dirname "$0")/../external/aocommon/scripts/format.sh
