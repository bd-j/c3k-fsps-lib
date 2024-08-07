#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Use the full resolution c3k library to produce smoothed interpolated on to
the logt-logg grid for fsps.
"""

import sys, os, shutil
from itertools import product
import logging

import numpy as np
import h5py
from scipy import interpolate

from prospect.sources import StarBasis
from utils_resample import make_seds
from utils_fsps import get_kiel_grid, interpolate_to_grid
from utils_fsps import get_binary_spec, interpolate_to_basel
from utils import get_ckc_parser, sed_to_bin
import constants


__all__ = ["sed", "rectify_sed", "to_grid", "to_bin", "make_fsps_metadata"]


logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger("resampler")
logger.setLevel(logging.INFO)


template = "{}/{}/{}_feh{:+3.2f}_afe{:+2.1f}.{}.h5"


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
    Generates an HDF5 file at `outname` continaing all resampled spectra for this afe-feh pair.

    outname : string
        The name of the HDF5 file into which the resampled spectra have been stored
    """
    specname = template.format(args.fulldir, args.specdir, args.ck_vers, feh, afe, "spec")
    fluxname = template.format(args.fulldir, args.fluxdir, args.ck_vers, feh, afe, "flux")
    outname = template.format(args.seddir, "", args.ck_vers, feh, afe, args.sedname)

    logger.info(f"reading spec from {specname}")
    logger.info(f"reading flux from {fluxname}")
    logger.info(f"putting output at {outname}")

    # Read Files and make the sed file
    specfile = h5py.File(specname, "r")
    fluxfile = h5py.File(fluxname, "r")
    make_seds(specfile, fluxfile, outname=outname,
              segments=segments, oversample=args.oversample,
              logger=logger)
    return outname


def rectify_sed(sedfile, rectified):
    """Do all the crazy magic to get models that are better for interpolating
    to BaSeL params.  Necessary?
    """
    with h5py.File(sedfile, "r") as sed:
        assert "extended" not in sed
        params = np.array(sed["parameters"])
        spec = np.array(sed["spectra"])
        wave = np.array(sed["wavelengths"])
        res = np.array(sed["resolution"])
        attrs = dict(sed.attrs)

    tgrid = np.sort(np.unique(params["logt"]))
    ggrid = np.sort(np.unique(params["logg"]))

    fudged = np.zeros(len(params))

    # --- fudge some spectra ---
    new_pars, new_spec, new_fudge = [], [], []

    # at each logt copy lowest logg to the two lower loggs
    for logt in tgrid:
        this = params["logt"] == logt
        gmin = np.min(params[this]["logg"])
        ind = (params["logt"] == logt) & (params["logg"] == gmin)
        for newg in [gmin - 0.5, gmin-1.0]:
            if newg >= ggrid.min():
                newp = np.array(params[ind])
                newp["logg"] = newg
                news = spec[ind]
                new_pars.append(newp)
                new_spec.append(news)
                new_fudge.append(1)

    lowt = tgrid[:16]

    # at 3.8 < logt < 4 copy max(logg) (=5.0 hopefully) to logg<=5.5 if needed
    for logt in tgrid:
        if (logt > 4.0) or (logt <= lowt.max()):
            continue
        this = params["logt"] == logt
        gmax = np.max(params[this]["logg"])
        ind = (params["logt"] == logt) & (params["logg"] == gmin)
        gtarg = np.arange(gmax+0.5, 5.51, 0.5)
        for newg in gtarg:
            newp = np.array(params[ind])
            newp["logg"] = newg
            news = spec[ind]
            new_pars.append(newp)
            new_spec.append(news)
            new_fudge.extend([1])

    # at all logt, interpolate across holes bracketed in logg
    for logt in tgrid:
        this = params["logt"] == logt
        x = params[this]["logg"]
        y = spec[this, :]
        newg = [logg for logg in ggrid if
                (logg > x.min()) & (logg < x.max()) & (logg not in x)]
        if len(newg) == 0:
            continue
        newg = np.array(newg)
        f = interpolate.interp1d(x, y, axis=0, kind="linear")
        news = f(newg)
        newp = np.tile(params[this][0], len(newg))
        newp["logg"] = newg
        new_pars.append(newp)
        new_spec.append(news)
        new_fudge.extend(len(newg) * [2])

    # at low logt (<~3.8), NN extrapolate from lower logg if missing
    for logt in lowt:
        this = params["logt"] == logt
        x = params[this]["logg"]
        y = spec[this, :]
        newg = [logg for logg in ggrid if
                (logg > x.max()) & (logg not in x)]
        if len(newg) == 0:
            continue
        newg = np.array(newg)
        f = interpolate.interp1d(x, y, axis=0, kind="nearest", fill_value="extrapolate")
        news = f(newg)
        newp = np.tile(params[this][0], len(newg))
        newp["logg"] = newg
        new_pars.append(newp)
        new_spec.append(news)
        new_fudge.extend(len(newg) * [1])

    if False:
        # at logt > 3.8, interpolate across holes in logg
        hight = tgrid[16:]
        for logg in ggrid:
            this = params["logg"] == logg
            x = params[this]["logt"]
            newt = [logt for i, logt in enumerate(hight[1:-1]) if
                    (hight[i-1] in x) & (hight[i+1] in x) and (logt not in x)]
            new_fudge.extend(len(newt) * [3])

    new_pars = np.concatenate(new_pars)
    new_spec = np.concatenate(new_spec)
    new_fudge = np.array(new_fudge)

    # concatenate the new models to the old ones
    params = np.concatenate([params, new_pars])
    spec = np.concatenate([spec, new_spec])
    fudged = np.concatenate([fudged, new_fudge])

    # write output to a file
    with h5py.File(rectified, "x") as sed:
        sed.create_dataset("extended", data=fudged)
        if "wavelengths" not in sed:
            sed.create_dataset("wavelengths", data=wave)
        if "resolution" not in sed:
            sed.create_dataset("resolution", data=res)
        #if "parameters" in sed:
        #    del sed["parameters"]
        sed.create_dataset("parameters", data=params)
        #if "spectra" in sed:
        #    del sed["spectra"]
        sed.create_dataset("spectra", data=spec)
        sed.attrs.update(attrs)


def to_grid(feh, afe, args, sedfile=None):
    """Interpolate spectrally resampled spectra for a given feh-afe pair onto
    the logt-logg grid and write to a new HDF5 file with extension `.fsps.h5`

    Parameters
    ----------
    feh : float
        Value of [Fe/H]

    afe : float
        Value of [alpha/Fe]

    args : namespace
        Must have the attributes
        * seddir - output location
        * sedname - add this label to the output filenames
        * use_basel_grid - whether to use the basel logt-logg grid or the one
                           given in ../data/log*.dat
        * nowrite - whether to write the output to HDF5 or return the
                    interpolated spectra

    sedfile : string
        Optional. Name of the HDF5 file with the resampled spectra for this
        feh-afe pair, if not given assumed to be
        `{seddir}/{ck_vers}_feh{feh:+3.2f}_afe{afe:+2.1f}.{sedname}.h5`
    """
    # Filenames
    sname = template.format(args.seddir, "", args.ck_vers, feh, afe, args.sedname)
    if sedfile is None:
        sedfile = sname
    outname = sedfile.replace(".h5", ".fsps.h5")

    # Grid Params and valid z=0.0200 spectra
    grid_pars = get_kiel_grid(basel=args.use_basel_grid)
    if args.use_basel_grid:
        _, _, valid = get_binary_spec(len(grid_pars), zstr="0.0200",
                                      speclib='BaSeL3.1/basel')
        interp_meth = interpolate_to_basel
    else:
        valid = np.ones(len(grid_pars), dtype=bool)
        interp_meth = interpolate_to_grid

    # Rectify?
    if getattr(args, "rectify", False):
        logger.info("Dilating C3K grid in logg before interpolating to FSPS grid")
        rectified = sedfile.replace(".h5", ".rect.h5")
        rectify_sed(sedfile, rectified)
        sedfile = rectified

    # Do the interpolation
    interpolator = StarBasis(sedfile, use_params=['logg', 'logt'], logify_Z=False,
                             n_neighbors=1, verbose=args.verbose,
                             rescale_libparams=True)
    bwave, bspec, inds = interp_meth(grid_pars, interpolator, valid=valid)

    # Keep track of how interpolation was done
    outside, inside, extreme = inds
    exact = (valid & (~outside) & (~inside) & (~extreme))

    if args.nowrite:
        return grid_pars, bwave, bspec, inds
    # Write the output
    logger.info(f"Writing interpolated spectra to {outname}")
    with h5py.File(outname, "w") as f:
        f.create_dataset("parameters", data=grid_pars)
        f.create_dataset("spectra", data=bspec)
        f.create_dataset("wavelengths", data=interpolator.wavelengths)
        idat = f.create_group("interpolation_info")
        idat.create_dataset("interpolated", data=inside)
        idat.create_dataset("exact", data=exact)
        idat.create_dataset("extrapolated", data=extreme)


def to_bin(feh=0.0, afe=0.0, zsol=0.0134,
           args=None, bindir=None, prefix=None):

    if bindir:
        args.bindir = bindir
    if prefix:
        args.prefix = prefix

    if getattr(args, "oldz", False):
        z = zsol * 10**feh
        zstr = f"{z:1.4f}"
        binname = f"{args.bindir}/{args.prefix}_z{z:1.4f}.spectra.bin"
    else:
        binname = f"{args.bindir}/{args.prefix}_feh{feh:+3.2f}_afe{afe:+2.1f}.spec.bin"

    sedfile = template.format(args.seddir, "", args.ck_vers, feh, afe, args.sedname)
    sedfile = sedfile.replace(".h5", ".fsps.h5")
    assert os.path.exists(sedfile)
    sed_to_bin(sedfile, binname)


def make_fsps_metadata(fehlist, args, zsol=0.0134):
    lname = f"{args.bindir}/{args.prefix}.lambda"
    zname = f"{args.bindir}/{args.prefix}_zlegend.dat"
    rname = f"{args.bindir}/{args.prefix}.res"
    wname = f"{args.bindir}/{args.prefix}.wave"

    # make the zlegend file
    zname = f"{args.bindir}/{args.prefix}_zlegend.dat"
    zlegend = open(zname, "w")
    for feh in np.sort(fehlist):
        if getattr(args, "oldz", False):
            z = zsol * 10**feh
            zstr = f"{z:1.4f}"
        else:
            zstr = f"{feh:+3.2f}"
        zlegend.write("{}\n".format(zstr))
    zlegend.close()

    # Now make the wavelength file
    sedfile = template.format(args.seddir, "", args.ck_vers, 0.0, 0.0, args.sedname)
    with h5py.File(sedfile, "r") as f:
        wave = np.array(f["wavelengths"])
        res = f["resolution"][:]
        vres = constants.ckms / res / constants.sigma_to_fwhm
        indef = (wave < 1e3) | (wave > 2.5e4)
        vres[indef] = -1 * vres[indef]
    wfile = open(wname, "w")
    lfile = open(lname, "w")
    rfile = open(rname, "w")
    for w, vr in zip(wave, vres):
        wfile.write("{}\n".format(w))
        rfile.write("{}\n".format(vr))
        lfile.write("{}  {}\n".format(w, vr))
    wfile.close()
    lfile.close()
    rfile.close()

    # now make sps_vars.f90
    return


if __name__ == "__main__":

    # These are the set of feh and afe from which we will choose based on zindex
    # These can be *OVERIDDEN* in a segments file.
    fehlist = [-2.0, -1.75, -1.5, -1.25, -1.0,
               -0.75, -0.5, -0.25, 0.0, 0.25, 0.5]
    afelist = [-0.2, 0.0, +0.2, +0.4, +0.6]

    parser = get_ckc_parser()
    parser.add_argument("--use_basel_grid", type=int, default=0)
    parser.add_argument("--rectify", type=int, default=1)
    parser.add_argument("--oldz", type=int, default=0)
    parser.add_argument("--make_seds", type=int, default=1)
    parser.add_argument("--make_grid", type=int, default=1)
    parser.add_argument("--make_bins", type=int, default=1)
    parser.add_argument("--force_meta", type=int, default=-1)
    parser.add_argument("--bindir", type=str, default="",
                        help="location of binary files")
    parser.add_argument("--nowrite", type=int, default=0)
    args = parser.parse_args()

    # -- Mess with some args ---
    args.fulldir = args.fulldir.format(args.ck_vers)
    logger.info(f"reading spectra from {args.fulldir}")
    if args.bindir == "":
        args.bindir = f"{args.seddir}/for_fsps"

    # --- read and copy segment file if it exists --
    if args.segment_file:
        logger.info(f"Reading segments from {args.segment_file}")
        import yaml
        import shutil
        with open(args.segment_file) as f:
            config = yaml.load(f, Loader=yaml.Loader)
        segments = config["segments"]
        if "fehlist" in config:
            fehlist = np.concatenate(config["fehlist"])
        if "afelist" in config:
            afelist = np.concatenate(config["afelist"])
        if "oversample" in config:
            args.oversample = config["oversample"]
            logger.info(f"Using oversample={args.oversample} from {args.segment_file}")
        else:
            logger.info(f"Using oversample={args.oversample}")
        if ("tag" in config) and (args.prefix == ""):
            args.prefix = f"{config['tag']}"

        _ = shutil.copy(args.segment_file, args.seddir)
        try:
            import json
            with open(os.path.join(args.segment_file, "args.json"), "w") as out:
                json.dump(args, out)
        except:
            pass

    # --- CHOOSE THE METALLICITY ---
    if args.zindex < 0:
        # for testing off odyssey
        feh, afe = 0.0, 0.0
        fehlist = [feh]
    else:
        metlist = list(product(fehlist, afelist))
        feh, afe = metlist[args.zindex]
    logger.info(f"Working on [Fe/H]={feh}, [alpha/Fe]={afe}")

    # --- make the sed file ---
    if args.make_seds:
        logger.info("Making H5 files with SEDs smoothed and stitched from native C3K spectral resolution")
        sedfile = sed(feh, afe, segments, args)

    # --- Make the SED interpolated to FSPS logt, logg grid ---
    if args.make_grid:
        logger.info("Interpolating smoothed SEDs to the FSPS Kiel grid")
        out = to_grid(feh, afe, args)
        if args.nowrite:
            grid_pars, bwave, bspec, inds = out

    # --- Make the binary file from the interpolated grid ---
    major = args.ck_vers.split("v")[-1][0]
    if major == "1":
        zsol = 0.0134
    elif major == "2":
        zsol = 0.019
    if args.make_bins:
        logger.info(f"Making spectral binary file at {args.bindir} with prefix {args.prefix}")
        to_bin(feh=feh, afe=afe, args=args, zsol=zsol)
        if ((afe == 0) & (feh == 0)) | (args.force_meta == args.zindex):
            logger.info("Making FSPS metadata files")
            make_fsps_metadata(fehlist, args, zsol=zsol)