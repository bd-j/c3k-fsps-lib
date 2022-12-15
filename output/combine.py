#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import numpy as np
import matplotlib.pyplot as pl

import h5py

if __name__ == "__main__":
    files = glob.glob("c3k_R10k/*feh?[0-3]*_afe+0.0*h5")
    files.sort()
    params, spectra = [], []
    for f in files:
        with h5py.File(f, "r") as dat:
            params.append(dat["parameters"][:])
            spectra.append(dat["spectra"][:])
            wave = dat["wavelengths"][:]
            attrs = dict(dat.attrs.items())

    parameters = np.concatenate(params)
    spectra = np.concatenate(spectra)
    with h5py.File("c3k_v1.3_afe+0.0_R10k.h5", "a") as out:
        out.create_dataset("parameters", data=parameters)
        out.create_dataset("spectra", data=spectra)
        out.create_dataset("wavelengths", data=wave)
        for k, v in attrs.items():
            out.attrs[k] = v