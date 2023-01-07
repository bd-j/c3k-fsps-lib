#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Use the full resolution c3k library to produce smoothed interpolated on to
the logt-logg grid for fsps.
"""

import sys, shutil
from itertools import product

import numpy as np
import h5py

from prospect.sources import StarBasis
from utils_resample import make_seds
from utils_fsps import get_kiel_grid, interpolate_to_grid
from utils import get_ckc_parser
    # key arguments are:
    #  * --zindex
    #  * --oversample
    #  * --ck_vers
    #  * --basedir
    #  * --seddir
    #  * --sedname

__all__ = ["sed", "to_grid"]


def sed(feh, afe, segments, args):
    """Resample the full resolution spectra for a single feh-afe pair to lower
    resolution spectra, and write to HDF5 file

    Parameters
    ----------
    feh : float
        Value of [Fe/H]

    afe : float
        Value of [alpha/Fe]

    segements : list of tuples
        Each tuple gives (wave_lo(float), wave_hi(float), R(float), FFT(bool))

    args : namespace
        Must have the attributes
        * ck_vers - version of C3K/CKC
        * fulldir - location of full resolution spectra HDF5 files
        * seddir - output location
        * sedname - add this label to the filenames
        * verbose
        * oversample - number of pixels per standard deviation of the LSF

    Returns
    -------
    Generates an HDF5 file at `outname`

    outname : string
        The name of the HDF5 file into which the resampled spectra have been stored
    """
    template = "{}/{}_feh{:+3.2f}_afe{:+2.1f}.{}.h5"
    specname = template.format(args.fulldir, args.ck_vers, feh, afe, "full")
    fluxname = template.format(args.fulldir, args.ck_vers, feh, afe, "flux")
    outname = template.format(args.seddir, args.ck_vers, feh, afe, args.sedname)

    # Read Files and make the sed file
    specfile = h5py.File(specname, "r")
    fluxfile = h5py.File(fluxname, "r")
    make_seds(specfile, fluxfile, fluxres=5e3, segments=segments,
              outname=outname, verbose=args.verbose,
              oversample=args.oversample)
    return outname


def to_grid(feh, afe, sedfile, args):
    """Interpolate spectrally resampled spectra onto the logt-logg grid,
    and write to a new HDF5 file with extension `.fsps.h5`

    Parameters
    ----------
    feh : float
        Value of [Fe/H]

    afe : float
        Value of [alpha/Fe]

    sedfile : string
        Name of the HDF5 file with the resampled spectra for this feh-afe pair

    args : namespace
        Must have the attributes
        * seddir - output location
        * sedname - add this label to the output filenames
        * nowrite - whether to write the output to HDF5 or return the
                    interpolated spectra
    """
    # Filenames
    template = "{}/{}_feh{:+3.2f}_afe{:+2.1f}.{}.fsps.h5"
    outname = template.format(args.seddir, args.ck_vers, feh, afe, args.sedname)

    # Grid Params and valid z=0.0200 spectra
    grid_pars = get_kiel_grid()
    valid = np.ones(len(grid_pars), dtype=bool)
    #cwave, cspec, valid = get_binary_spec(len(grid_pars), zstr="0.0200",
    #                                      speclib='BaSeL3.1/basel')

    # My interpolator
    interpolator = StarBasis(sedfile, use_params=['logg', 'logt'], logify_Z=False,
                             n_neighbors=1, verbose=args.verbose,
                             rescale_libparams=True)

    # Do the interpolation
    bwave, bspec, inds = interpolate_to_grid([grid_pars, cwave, cspec], interpolator,
                                             valid=valid, renorm=False, plot=None,
                                             verbose=args.verbose)
    # Keep track of how interpolation was done
    false = np.zeros(len(grid_pars), dtype=bool)
    o, i, e = inds
    out, interp, extreme = false.copy(), false.copy(), false.copy()
    out[o] = True
    interp[i] = True
    extreme[e] = True
    exact = (valid & (~out) & (~interp) & (~extreme))

    if args.nowrite:
        return grid_pars, bwave, bspec, inds
    # Write the output
    with h5py.File(outname, "w") as f:
        f.create_dataset("parameters", data=grid_pars)
        f.create_dataset("spectra", data=bspec)
        f.create_dataset("wavelengths", data=interpolator.wavelengths)
        idat = f.create_group("interpolation_info")
        idat.create_dataset("interpolated", data=interp)
        idat.create_dataset("exact", data=exact)
        idat.create_dataset("nearest_tg", data=extreme)


if __name__ == "__main__":

    # These are the set of feh and afe from which we will choose based on zindex
    # These can be *OVERIDDEN* in a segments file.
    fehlist = [-2.0, -1.75, -1.5, -1.25, -1.0,
               -0.75, -0.5, -0.25, 0.0, 0.25, 0.5]
    afelist = [-0.2, 0.0, +0.2, +0.4, +0.6]

    parser = get_ckc_parser()
    parser.add_argument("--segment_file", type=str, default=None)
    parser.add_argument("--nowrite", type=int, default=0)
    args = parser.parse_args()

    # -- Mess with some args ---
    args.fulldir = args.fulldir.format(args.ck_vers)

    # --- read and copy segment file if it exists --
    import yaml
    import shutil
    with open(args.segment_file) as f:
        config = yaml.load(f, Loader=yaml.Loader)
    segments = config["segments"]
    try:
        fehlist = np.concatenate(config["fehlist"])
    except(KeyError):
        print("Using default fehlist")
    try:
        afelist = np.concatenate(config["afelist"])
    except(KeyError):
        print("Using default afelist")
    if "oversample" in config:
        args.oversample = config["oversample"]
        print(f"Using oversample={args.oversample} from {args.segment_file}")
    else:
        print(f"Using oversample={args.oversample}")

    _ = shutil.copy(args.segment_file, args.seddir)

    # --- CHOOSE THE METALLICITY ---
    if args.zindex < 0:
        # for testing off odyssey
        feh, afe = 0.0, 0.0
    else:
        metlist = list(product(fehlist, afelist))
        feh, afe = metlist[args.zindex]
    print(feh, afe)

    # --- make the sed file ---
    sedfile = sed(feh, afe, segments, args)

    # --- Make the SED interpolated to FSPS logt, logg grid ---
    if args.nowrite:
        sys.exit()
    if "sedfile" not in locals():
        template = "{}/{}_feh{:+3.2f}_afe{:+2.1f}.{}.h5"
        sedfile = template.format(args.seddir, args.ck_vers, feh, afe, args.sedname)
    out = to_grid(feh, afe, sedfile, args)
    if args.nowrite:
        grid_pars, bwave, bspec, inds = out
