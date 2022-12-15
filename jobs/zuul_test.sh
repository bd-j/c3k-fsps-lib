#!/usr/bin/env bash

libname=nirspec
python c3k_resample.py --segment_file ../segments/segments_${libname}.yml --oversample 2 \
                       --seddir ../output/${libname} --sedname ${libname} \
                       --fulldir "/Users/bjohnson/Projects/ckc/ckc/spectra/fullres/c3k/"