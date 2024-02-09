#!/bin/bash

# Optional command line parameters:
# $1: number of samples per subpixel (non-uniform filtering) - called "samples" in the program "PathTracing.cpp".
# $2: aperture size
# $3: focal length

rm -r build/
mkdir build/
cd build/
cmake -DCMAKE_BUILD_TYPE=RELEASE ../
make
cd bin/
./PathTracing $1 $2 $3 