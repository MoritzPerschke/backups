#!/bin/bash

# Optional command line parameters:
# $1: number of samples per subpixel (non-uniform filtering) - called "samples" in the program "PathTracing.cpp".
# $2: aperture size
# $3: focal length

echo '[*] Running Pathtracing executable:'
cd build/
make
cd bin/
./PathTracing $1 $2 $3 
echo '[+] Render done, exiting...'
