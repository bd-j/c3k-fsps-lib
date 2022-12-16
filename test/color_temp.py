
import argparse
import matplotlib.pyplot as pl
import numpy as np
import h5py
from sedpy import observate

caa = 2.998e18

def get_mags(fsps_h5, filternames):
    specfile = fsps_h5
    wave = specfile['wavelengths'][:]
    spec = specfile["spectra"][:]
    spec *= caa / wave**2
    filters = observate.FilterSet(filternames)
    mags = observate.getSED(wave, spec, filters)

    return mags, specfile["parameters"][:]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--specfile", type=str, default="../output/nirspec/c3k_v1.3_feh+0.00_afe+0.0.nirspec.fsps.h5")
    args = parser.parse_args()

    filternames = [f"sdss_{b}0" for b in "ugriz"]
    specfile = h5py.File(args.specfile, "r")
    mags, params = get_mags(specfile, filternames)
    sel = specfile["spectra"][:, 3000] > 1e-33

    pl.ion()
    fig, ax = pl.subplots()
    cb = ax.scatter(params["logt"][sel], mags[sel, 1] - mags[sel, 3], c=params["logg"][sel])
    ax.set_xlabel(r"$\log {\rm T}_{eff}$")
    ax.set_ylabel(r"$g - i$")
    ax.text
    pl.colorbar(cb, label=r"$\log g$")