#!/usr/bin/env bash

libname=lr
mkdir -p ../output/${libname}
mkdir -p ../output/${libname}/for_fsps

python c3k_resample.py --segment_file ../segments/segments_${libname}.yml --zindex -99 --oversample 2 \
                       --seddir ../output/${libname} --sedname ${libname} \
                       --fulldir "/Users/bjohnson/Projects/ckc/ckc/spectra/fullres/c3k/"