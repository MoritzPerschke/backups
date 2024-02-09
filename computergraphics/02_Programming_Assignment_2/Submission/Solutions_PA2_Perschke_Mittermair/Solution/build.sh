#!/bin/bash

echo '[*] Building Pathtracing executable'
rm -r build/
mkdir build/
cd build/
cmake -DCMAKE_BUILD_TYPE=RELEASE ../
make  
echo '[+] Done...'
