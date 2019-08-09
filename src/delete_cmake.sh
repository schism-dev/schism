#!/bin/sh
# Eliminate cmake files in the /src directory.

rm -f cmake_install.cmake CMakeCache.txt
find . -name CMakeFiles -type d -print | xargs rm -rf
