#!/bin/bash

cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ ..
make
rm -r CMakeFiles
