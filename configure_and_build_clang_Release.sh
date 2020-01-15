#!/bin/bash
# Please do not modify this script

# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact with us in:
#   - https://github.com/KratosMultiphysics/Kratos

# Optional parameters:
# You can find a list will all the compiation options in INSTALL.md or here:
#   - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
#export CC=/usr/bin/clang-9
#export CXX=/usr/bin/clang++-9

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:-"Release"}

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/ExternalSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/MultiscaleROMApplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
-D CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 " \
-D CMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 " \
-D USE_MPI=OFF \
-D KRATOS_EXCLUDE_OPENMP=ON \
-D INSTALL_RUNKRATOS=OFF \
-D KRATOS_BUILD_TESTING=ON \
-D KRATOS_COLORED_OUTPUT=OFF \

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j 7
