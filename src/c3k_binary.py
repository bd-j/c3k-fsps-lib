#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Method and for converting individual metallicity (feh/afe combo) .fsps.h5
HDF5 files of SEDs to FSPS compatible binary format with ancillary info files.

This assumes the SEDs have already been interpolated to the BaSeL stellar
parameter grid, and that the spectra in the HDF5 are ordered with logt changing
fastest and then logg.
"""

import glob
import numpy as np
import matplotlib.pyplot as pl
import h5py
from argparse import Namespace

from utils import get_ckc_parser, sed_to_bin
import constants

__all__ = ["prep_for_fsps"]


def prep_for_fsps(fehlist=[-0.75, -0.5, -0.25, 0.0, 0.25],
                  zsol=0.0134, afe=0.0, args=None, outdir=None, prefix=None):

    if outdir is not None:
        args.outdir = outdir
    if prefix is not None:
        args.prefix = prefix

    template = "{}/{}_feh{:+3.2f}_afe{:+2.1f}.{}.fsps.h5"
    binary_name = "{}/{}_z{}.spectra.bin"
    lambda_name = "{}/{}.lambda"
    zname = "{}/{}_zlegend.dat"
    rname = "{}/{}.res"

    zlegend = open(zname.format(args.outdir, args.prefix), "w")
    for j, feh in enumerate(fehlist):
        z = zsol * 10**feh
        zstr = "{:1.4f}".format(z)
        zlegend.write("{}\n".format(zstr))
        sedfile = template.format(args.seddir, args.ck_vers, feh, afe, args.sedname)
        outname = binary_name.format(args.outdir, args.prefix, zstr)
        sed_to_bin(sedfile, outname)
    zlegend.close()
    # Now make the wavelength file
    with h5py.File(sedfile.replace(".fsps.", "."), "r") as f:
        wave = np.array(f["wavelengths"])
        res = f["resolution"][:]
        vres = constants.ckms / res / constants.sigma_to_fwhm
        indef = (wave < 1e3) | (wave > 2.5e4)
        vres[indef] = -1 * vres[indef]
    with open(lambda_name.format(args.outdir, args.prefix), "w") as wavefile:
        for w in wave:
            wavefile.write("{}\n".format(w))
    with open(rname.format(args.outdir, args.prefix), "w") as resfile:
        for vr in vres:
            resfile.write("{}\n".format(vr))


if __name__ == "__main__":

    parser = get_ckc_parser()
    parser.add_argument("--outdir", type=str, default="")
    parser.add_argument("--test", type=int, default=0)

    args = parser.parse_args()

    afelist = [-0.2, 0.0, 0.2, 0.4, 0.6]
    fehlist = [-2.0, -1.75, -1.5, -1.25, -1.0,
               -0.75, -0.5, -0.25, 0.0, 0.25, 0.5]

    afe = afelist[args.zindex]
    if args.test:
        fehlist = [0.0]

    # should start with c3k or sps_setup.f90 needs to change
    args.prefix = "c3k_afe{:+2.1f}".format(afe)

    prep_for_fsps(fehlist=fehlist, afe=afe, args=args)
