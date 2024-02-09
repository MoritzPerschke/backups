#!/bin/bash

# Optional command line parameter: number of samples per subpixel (non-uniform filtering) - called "samples" in the program "PathTracing.cpp".

cd build/
make
cd bin/
./PathTracing $1