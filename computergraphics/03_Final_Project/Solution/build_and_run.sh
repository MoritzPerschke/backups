#!/bin/bash

# Optional command line parameter: number of samples per subpixel (non-uniform filtering) - called "samples" in the program "PathTracing.cpp".

rm -r build/
mkdir build/
cd build/
cmake -DCMAKE_BUILD_TYPE=RELEASE ../
make
cd bin/
./PathTracing $1