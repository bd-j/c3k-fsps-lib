#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module for converting single metallicity SED files in HDF format to the
BaSeL parameter grid expected by fsps.  Lots of interpolation.
"""

import os
from itertools import product
import numpy as np
import h5py

from scipy import interpolate

from utils import read_binary_spec
from prospect.sources import StarBasis

__all__ = ["get_kiel_grid", "get_binary_spec", "interpolate_to_grid"]


def dict_struct(strct):
    """Convert from a structured array to a dictionary.  This shouldn't really
    be necessary.
    """
    return dict([(n, strct[n]) for n in strct.dtype.names])


def get_kiel_grid(basel=False):
    """Get the logt-logg grid parameters as a 1-d list of parameter tuples.  The
    binary files are written in the order logg, logt, wave with wave changing
    fastest and logg the slowest.
    """
    if basel:
        dirn = os.path.join(os.environ["SPS_HOME"], "SPECTRA", "BaSeL3.1")
        pre = "basel_"
    else:
        dirn = "../data"
        pre = ""

    logg = np.genfromtxt(f"{dirn}/{pre}logg.dat")
    logt = np.genfromtxt(f"{dirn}/{pre}logt.dat")
    logt = np.log10(np.round(10**logt))
    ngrid = len(logg) * len(logt)
    dt = np.dtype([('logg', np.float), ('logt', np.float)])
    kiel_params = np.array(list(product(logg, logt)), dtype=dt)
    return kiel_params


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
    new_pars, new_spec = [], []

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

    # at lowest and second lowest logt, interpolate across holes in logg
    for logt in tgrid[:2]:
        this = params["logt"] == logt
        x = params[this]["logg"]
        y = spec[this, :]
        newg = [logg for logg in ggrid if
                (logg > x.min()) & (logg < x.max()) and (logg not in x)]
        if len(newg) == 0:
            continue

        newg = np.array(newg)
        f = interpolate.interp1d(x, y, axis=0, kind="linear")
        news = f(newg)
        newp = np.tile(params[this][0], len(newg))
        newp["logg"] = newg
        new_pars.append(newp)
        new_spec.append(news)

    new_pars = np.concatenate(new_pars)
    new_spec = np.concatenate(new_spec)

    # concatenate the new models to the old ones
    params = np.concatenate([params, new_pars])
    spec = np.concatenate([spec, new_spec])
    fudged = np.concatenate([fudged, np.ones(len(new_pars))])

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


def interpolate_to_grid(grid_pars, interpolator, valid=True, verbose=True):

    bwave = interpolator.wavelengths
    allspec = []
    outside, inside, extreme = [], [], []

    for i, (p, v) in enumerate(zip(grid_pars, valid)):
        inds, wghts = interpolator.weights(**dict_struct(p))
        ex_g = extremeg(interpolator, p)
        if verbose:
            print(i, p, v, ex_g, len(inds))
        if ex_g and v:
            inds, wghts = nearest_tg(interpolator, p)

        # valid *and* (in-hull and non-extreme g)
        if (v and (len(inds) > 1)):
            _, bspec, _ = interpolator.get_star_spectrum(**dict_struct(p))
        # valid *and* (out-hull or extremeg)
        elif (v and (len(inds) == 1)):
            bspec = interpolator._spectra[inds, :]
        # not valid
        else:
            bspec = np.zeros_like(interpolator.wavelengths) + 1e-33
        allspec.append(np.squeeze(bspec))

    return interpolator.wavelengths, np.array(allspec), [outside, inside, extreme]


def nearest_tg(interpolator, pars):
    """Do nearest neighbor but, first find the nearest logt and only then the
    nearest logg
    """
    tgrid = np.unique(interpolator._libparams["logt"])
    tt = tgrid[np.argmin(abs(tgrid - pars["logt"]))]
    thist = interpolator._libparams["logt"] == tt
    ggrid = np.unique(interpolator._libparams["logg"][thist])
    gg = ggrid[np.argmin(abs(ggrid - pars["logg"]))]
    choose = ((interpolator._libparams["logt"] == tt) &
              (interpolator._libparams["logg"] == gg))
    assert choose.sum() == 1
    ind = np.where(choose)[0][0]
    wght = 1.0
    return np.array([ind]), np.array([wght])


def extremeg(interpolator, pars):
    """Find if a set of parameters is below the lowest existing gravity for
    that temperature.
    """
    tgrid = np.unique(interpolator._libparams["logt"])
    tt = tgrid[np.argmin(abs(tgrid - pars["logt"]))]
    thist = interpolator._libparams["logt"] == tt
    ggrid = np.unique(interpolator._libparams["logg"][thist])
    return (pars["logg"] < ggrid.min()) or (pars["logg"] > ggrid.max())


def get_binary_spec(ngrid, zstr="0.0200", speclib="BaSeL3.1/basel"):
    """
    :param zstr: for basel "0.0002", "0.0006", "0.0020", "0.0063", "0.0200", "0.0632"
    """
    specname = "{}/SPECTRA/{}".format(os.environ["SPS_HOME"], speclib)
    wave = np.genfromtxt("{}.lambda".format(specname))
    try:
        ss = read_binary_spec("{}_wlbc_z{}.spectra.bin".format(specname, zstr), len(wave), ngrid)
    except(IOError):
        ss = read_binary_spec("{}_z{}.spectra.bin".format(specname, zstr), len(wave), ngrid)
    #logg = np.genfromtxt("{}_logg.dat".format(specname))
    #logt = np.genfromtxt("{}_logt.dat".format(specname))
    #spec = ss.reshape(len(logg), len(logt), len(wave))
    valid = ss.max(axis=1) > 1e-32
    return wave, ss, valid


def renorm(spec, normed_spec, wlo=5e3, whi=2e4):

    w, f = spec
    wn, fn = normed_spec
    f = np.squeeze(f)
    fn = np.squeeze(fn)
    g = (w > wlo) & (w < whi)
    l = np.trapz(f[g], w[g])
    g = (wn > wlo) & (wn < whi)
    ln = np.trapz(fn[g], wn[g])
    return ln / l, f * ln / l


def comp_text(inds, wghts, interpolator):
    txt = ""
    for i, w in zip(inds, wghts):
        cd = dict_struct(interpolator._libparams[i])
        txt += "\n{w:0.2f}@{i}: {logt:4.3f} {logg:3.2f}".format(i=i, w=w, **cd)

    return txt


def plot_interp(cwave, cflux, bflux, inds, wghts, interpolator,
                show_components=True, renorm_pars={}):

    import matplotlib.pyplot as pl
    bwave = interpolator.wavelengths
    _, bs = renorm([bwave, bflux], [cwave, cflux], **renorm_pars)
    fig, ax = pl.subplots()
    ax.plot(cwave[::5], cflux[::5], label="Charlie binary")
    ax.plot(bwave, bs, label="interpolated")
    ax.set_xlim(1e3, 3e4)

    if show_components:
        txt = comp_text(inds, wghts, interpolator)
        for i, w in zip(inds, wghts):
            if w > 0:
                _, bs = renorm([bwave, interpolator._spectra[i, :]],
                                [cwave, cflux], **renorm_pars)
                ax.plot(bwave, bs, label="comp. {}, wght={}".format(i, w), alpha=0.5)
        ax.text(0.3, 0.3, txt, transform=ax.transAxes)

    return fig, ax


def compare_at(charlie, interpolator, logg=4.5, logt=np.log10(5750.),
               show_components=False, renorm_pars={}):

    bpars, cwave, cspec = charlie
    sel = (bpars["logg"] == logg) & (bpars["logt"] == logt)
    if sel.sum() != 1:
        print("This set of parameters is not in the BaSeL grid: logg={}, logt={}".format(logg, logt))
        raise(ValueError)

    cflux = np.squeeze(cspec[sel, :])
    assert cflux.max() > 1e-33, ("requested spectrum has no "
                                 "Charlie spectrum: logg={}, logt={}".format(logg, logt))

    bwave = interpolator.wavelengths
    bflux, _, _ = interpolator.get_spectrum(logg=logg, logt=logt)
    inds, wghts = interpolator.weights(logg=logg, logt=logt)

    fig, ax = plot_interp(cwave, cflux, bflux, inds, wghts, interpolator)
    ax.set_title("target: {logt:4.3f} {logg:3.2f}".format(logg=logg, logt=logt))

    return fig, ax


def show_coverage(grid_pars, libparams, inds, valid):
    false = np.zeros(len(grid_pars), dtype=bool)
    o, i, e = inds
    out, interp, extreme = false.copy(), false.copy(), false.copy()
    out[o] = True
    interp[i] = True
    extreme[e] = True
    exact = (valid & (~out) & (~interp) & (~extreme))

    import matplotlib.pyplot as pl
    fig, ax = pl.subplots()
    ax.plot(grid_pars["logt"], grid_pars["logg"], 'o', alpha=0.2,
            label="BaSeL grid")
    ax.plot(grid_pars["logt"][exact], grid_pars["logg"][exact], 'o',
            label="exact")
    ax.plot(grid_pars["logt"][out], grid_pars["logg"][out], 'o',
            label="outside C3K")
    ax.plot(grid_pars["logt"][interp], grid_pars["logg"][interp], 'o',
            label="interpolated")
    ax.plot(grid_pars["logt"][extreme], grid_pars["logg"][extreme], 'o',
            label="nearest (t, g)")
    ax.plot(libparams["logt"], libparams["logg"], 'ko', alpha=0.3, label="C3K")
    ax.invert_yaxis()
    ax.invert_xaxis()
    #ax.legend(loc=0)

    return fig, ax
