#!/bin/bash
cd build
#cmake -DCMAKE_BUILD_TYPE=Debug ..
#cmake -DCMAKE_CXX_FLAGS="-DOPT_HIT_ONCE" ..
#cmake -DCMAKE_CXX_FLAGS="-DOPT_CHECK_ONCE" ..
cmake ..
make -j
