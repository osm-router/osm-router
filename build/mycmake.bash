#!/bin/bash

cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++-3.5 ..
make
rm -r CMakeFiles
