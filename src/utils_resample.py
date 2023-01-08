# -*- coding: utf-8 -*-


"""Module for producing full range spectra from c3k .spec and .flux files (or
rather their hdf5 versions) This involves convolving the appropriate spectra by
the appropriate amounts for a set of segements and stitching them together.
Note that the .flux files actually are top-hat filtered versions of the spectra,
with the value at each point giving the sum of total flux within some +/- from
that point, and have some weird effective resolution.
"""
import numpy as np
import h5py, json
from sedpy.smoothing import smoothspec
from utils import segments_to_wavelength, construct_outwave, constants

__all__ = ["make_seds", "make_one_sed"]


def make_seds(specfile, fluxfile, outname=None,
              segments=[()], oversample=2,
              specres=1e5, fluxres=5e3, logger=None):
    """
    Parameters
    -------
    specfile : Dict-like with keys "parameters", "wavelengths", and "spectra"
        Handle to the HDF5 file containting the high-resolution C3K spectra

    fluxfile : Dict-like with keys "parameters", "wavelengths", and "spectra"
        Handle to the HDF5 file containing the low-resolution C3K "flux"
        spectra, which is assumed matched line by line to the specfile

    outname : string (optional, default: None)
        Name of output h5 file.  If not given, then method will return
        parameters, interpolated spectra, and a tuple of wavelength and
        resolution.

    segments : list of 4-tuples
        A list of 4-tuples describing the wavelength segments and their
        resolution, and whether to use FFT (not recommended near the edge o the
        hires file).  Each tuple should have the form (lo, hi, R, bool)

    oversample : integer, (default: 2)
        Number of pixels per FWHM

    Returns
    ------
    Generates an HDF file at `outname`
    """
    if logger is None:
        log = print
    else:
        log = logger.info
    # --- Wavelength arrays ----
    swave = np.array(specfile["wavelengths"])
    fwave = np.array(fluxfile["wavelengths"])
    outwave, outres = segments_to_wavelength(segments, oversample=oversample)
    assert np.all(np.diff(outwave) > 0), "Output wavelength grid is not ascending!"
    nw = len(outwave)

    # --- Match specfile to fluxfile ----
    # this is probably a stupid way to do this.
    sind, find = [], []
    for i, spec in enumerate(specfile["spectra"]):
        # find the matching flux entry
        params = specfile["parameters"][i]
        ind = [fluxfile["parameters"][f] == params[f] for f in params.dtype.names]
        ind = np.array(ind).prod(axis=0).astype(int)
        if ind.sum() != 1:
            pdict = dict(zip(constants.param_order, params))
            log(f"could not find unique flux spectrum @ params {pdict}")
        else:
            sind.append(i)
            find.append(ind.tolist().index(1))
    sind = np.array(sind)
    find = np.array(find)

    # --- Setup the output h5 file ---
    if outname is not None:
        out = h5py.File(outname, "w")
        partype = specfile["parameters"].dtype
        wavout = out.create_dataset('wavelengths', data=outwave)
        resout = out.create_dataset('resolution', data=outres)
        sedout = out.create_dataset('spectra', shape=(0, nw), maxshape=(None, nw))
        parout = out.create_dataset('parameters', shape=(0,), maxshape=(None,), dtype=partype)
        out.attrs["segments"] = json.dumps(segments)
        out.attrs["segments_desc"] = "(lo, hi, R_fwhm, FFT)"
    else:
        nsed = len(sind)
        sedout = np.zeros([nsed, nw])
        parout = np.zeros([nsed], dtype=partype)

    #  --- Fill H5 file ---
    # loop over spectra convolving segments and getting the SEDs, and putting
    # them in the SED file
    #matches = zip(specfile["parameters"][sind],
    #              specfile["spectra"][sind, :],
    #              fluxfile["spectra"][find, :])
    for i, (s, f) in enumerate(zip(sind, find)):
        spec = specfile["spectra"][s, :]
        flux = fluxfile["spectra"][f, :]
        wave, sed, msg = make_one_sed(swave, spec, fwave, flux, segments,
                                      oversample=oversample,
                                      specres=specres, fluxres=fluxres)
        assert len(sed) == nw, (f"SED is not the same length as the desired "
                                f"output wavelength grid! ({len(sed)} != {len(outwave)})")

        sed[np.isnan(sed)] = 0.0
        if outname is not None:
            sedout.resize(i+1, axis=0)
            parout.resize(i+1, axis=0)
        sedout[i, :] = sed
        parout[i] = specfile["parameters"][s]
        if i == 0:
            log(msg)

    if outname:
        out.close()
        return
    else:
        return np.array(parout), np.array(sedout), (outwave, outres)


def make_one_sed(swave, spec, fwave, flux, segments=[()], clip=1e-33,
                 specres=1e5, fluxres=5000., oversample=2, verbose=False):
    """
    Parameters
    -----------
    swave : ndarray of shape (nwh,)
        Full resolution spectrum wavelength vector for the `.spec` output of synthe.

    spec : ndarray of shape (nws,)
        Full resolution flux vector.

    fwave : ndarray of shape (nwl,)
        Wavelength vector of the lower resolution `flux` spectrum provided by synthe.

    flux : ndarray of shape (nwl,)
        Low resolution flux vector, same units as `spec`

    segments: list of tuples
        Specification of the wavelength segments, of the form
        (wave_min, wave_max, resolution, use_fft)

    clip : float (default: 1e-33)
        Lower limit for the output SED fluxes

    specres : float (default: 3e5)
        The resolution of the input spectrum given by the `swave` and `spec` vectors.

    fluxres : float (default: 5e3)
        The resolution of the input spectrum given by the `fwave` and `flux` vectors.

    oversample : float (default, 2)
        Number of output pixels per FWHM of the line-spread function.

    Returns
    --------

    outwave : ndarray of shape (nout,)
        Wavelength vector of the output downsampled SED

    sed : ndarray of shape (nout,)
        Flux vector of the output downsampled SED (same units as `spec` and `flux`)
    """
    sed = []
    outwave = []
    for j, (lo, hi, rout, fftsmooth) in enumerate(segments):
        # get the output wavelength vector for this segment, throwing away
        # the last point (which will be the same as the first of the next
        # segment)
        out = construct_outwave(lo, hi, resolution=rout, logarithmic=True, oversample=oversample)[:-1]
        # Do we use the hires or lores spectrum?
        if (lo > swave.min()) and (hi < swave.max()):
            inspec = spec
            inres = specres
            inwave = swave
            msg = f"using hires for {lo} - {hi} @ R={rout}"
        else:
            inspec = flux
            inres = fluxres
            inwave = fwave
            msg = f"using lores for {lo} - {hi} @ R={rout}"

        if fftsmooth:
            msg += "; using FFT"
        if verbose:
            print(msg)
        assert rout <= inres, "You are trying to smooth to a higher resolution than than the input provides!"
        # account for C3K resolution
        rsmooth = (rout**(-2.) - inres**(-2))**(-0.5)
        # convert to lambda/sigma_lambda
        rsmooth *= constants.sigma_to_fwhm
        s = smoothspec(inwave, inspec, rsmooth, smoothtype="R", outwave=out, fftsmooth=fftsmooth)
        if clip > 0:
            np.clip(s, clip, np.inf, out=s)

        sed.append(s)
        outwave.append(out)

    outwave = np.concatenate(outwave)
    sed = np.concatenate(sed)
    # now replace the lambda > fwave.max() (plus some padding) with a BB
    # Assumes units are fnu
    fwave_max = np.max(fwave[flux > 0])
    ind_max = np.searchsorted(outwave, fwave_max)
    sed[ind_max-9:] = sed[ind_max - 10] * (outwave[ind_max - 10] / outwave[ind_max-9:])**2

    return outwave, sed, msg
