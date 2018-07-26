#!/bin/bash
cd isochrone
rm obj/*.o
rm lib/libisodist_js.so
rm isodist_js.so
make isodist_js.so -j 15
make lib/libisodist_js.so -j 15
cd ../edf/
rm lib/libedf.so
make lib/libedf.so -j 15
rm py/edf_py.so
make py/edf_py.so -j 15
cd ../edf_sampling/
rm py/edf_sampling.so
make py/edf_sampling.so
cd ..
