#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import struct
import h5py


__all__ = ["get_ckc_parser",
           "segments_to_wavelength", "construct_outwave",
           "convert_resolution",
           "read_binary_spec", "sed_to_bin"]

import constants


def segments_to_wavelength(segments, oversample=2):
    outwave = [construct_outwave(lo, hi, resolution=rout, logarithmic=True, oversample=oversample)[:-1]
               for (lo, hi, rout, _) in segments]
    res = [np.ones_like(w)*rout for w, (lo, hi, rout, _) in zip(outwave, segments)]
    return np.concatenate(outwave), np.concatenate(res)


def construct_outwave(min_wave_smooth=0.0, max_wave_smooth=np.inf,
                      dispersion=1.0, oversample=2.0,
                      resolution=3e5, logarithmic=False,
                      **extras):
    """Given parameters describing the output spectrum, generate a wavelength
    grid that properly samples the resolution.

    Parameters
    ----------
    min_wave_smooth : float, optional (default: 0.)
        Minimum wavelength of the output wavelength vector

    max_wave_smooth : float (default: np.inf)
        Maximum wavelength of the output wavelength vector,
        same units as min_wave_smooth

    dispersion : float, optional (default: 1)
        Wavelength units per resolution element if `logarithmic=False`

    resolution : float, optional
        The value of the wavelength divided by the resolution element,
        if `logarithmic=True`

    oversample : float, optional, (default: 2)
        The number of pixels per resolution element.

    Returns
    ---------
    wave : ndarray
        The output wavelength vector.
    """
    if logarithmic:
        # critically sample the resolution
        dlnlam = 1.0 / resolution / oversample
        lnmin, lnmax = np.log(min_wave_smooth), np.log(max_wave_smooth)
        #print(lnmin, lnmax, dlnlam, resolution, oversample)
        out = np.exp(np.arange(lnmin, lnmax + dlnlam, dlnlam))
    else:
        out = np.arange(min_wave_smooth, max_wave_smooth,
                        dispersion / oversample)
    return out


def convert_resolution(R_fwhm, R_library=3e5):
    """Convert from standard 'R" values based on lambda/FWHM to the
    lambda/sigma values expected by smoothspec
    """
    R_sigma = constants.sigma_to_fwhm / np.sqrt(1 / R_fwhm**2 - 1 / R_library**2)
    return R_sigma


def get_ckc_parser():
    """A general purpose argument parser for ckc methods and modules.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", type=bool, default=True,
                        help="chatter?")
    parser.add_argument("--np", type=int, default=1,
                        help="number of processors")

    # Metallicity choices
    parser.add_argument("--zindex", type=int, default=-99,
                        help="index of the metallicity mixture to use")
    parser.add_argument("--feh", type=float, default=0.0,
                        help=("The feh value to process."))
    parser.add_argument("--afe", type=float, default=0.0,
                        help=("The afe value to process."))

    # File specifying  resolution segments
    parser.add_argument("--segment_file", type=str, default=None)

    # Number of pixels per sigma
    parser.add_argument("--oversample", type=float, default=2,
                        help="Factor by which to oversample the dispersion")

    # Filenames, library versions, and spectrum locations
    parser.add_argument("--ck_vers", type=str, default="c3k_v2.2",
                        help=("Name of directory that contains the version of C3K spectra to use."))
    parser.add_argument("--basedir", type=str,
                        default='/n/holystore01/LABS/conroy_lab/Lab/cconroy/kurucz/grids',
                        help=("Location of the directories containing different C3K versions."))
    parser.add_argument("--fulldir", type=str,
                        default='/n/holystore01/LABS/conroy_lab/Lab/bdjohnson/data/kurucz/{}/fullopt/',
                        help=("Location to store the HDF5 versions of .spec and .flux"))
    parser.add_argument("--fluxdir", type=str, default="flux")
    parser.add_argument("--specdir", type=str, default="spec")
    parser.add_argument("--seddir", type=str, default='./',
                        help=("Path to the directory where the resampled sed files will be placed."))
    parser.add_argument("--sedname", type=str, default="lr",
                        help=("nickname for the SED file, e.g. sedR500"))
    parser.add_argument("--prefix", type=str, default="",
                        help=("prefix for the binary files"))

    # For `make_full_h5` if making HDF5 for the first time
    parser.add_argument("--spec_type", type=str, default='lores',
                        help=("Whether to create HDF5 for the full resolution "
                              "spectra ('hires') or for the low resolution, "
                              "larger wavelength range .flux files ('lores')"))

    return parser


def read_binary_spec(filename, nw, nspec):
    """Read a binary file with name ``filename`` containing ``nspec``
    spectra each of length ``nw`` wavelength points and return an
    array of shape (nspec, nw)
    """
    count = 0
    spec = np.empty([nspec, nw])
    with open(filename, 'rb') as f:
        while count < nspec:
                count += 1
                for iw in range(nw):
                    byte = f.read(4)
                    spec[count-1, iw] = struct.unpack('f', byte)[0]
    return spec


def sed_to_bin(sedfile, outname):
    """Note the sedfile is expected to already have the spectra in the correct
    order (logt changing fastest, then logg)
    """
    with h5py.File(sedfile, "r") as f:
        spectra = np.array(f["spectra"])
    with open(outname, "wb") as outfile:
        for spec in spectra:
            for flux in spec:
                if flux < 1e-33:
                    flux = 1e-33
                outfile.write(struct.pack('f', flux))


def combine_h5(files, outname):
    wave, parameters, spectra = None, [], []
    for f in files:
        with h5py.File(f, "r") as h5:
            parameters.append(h5["parameters"][:])
            spectra.append(h5["spectra"][:])
            w = h5["wavelengths"][:]
            if wave is None:
                wave = w
            else:
                assert np.allclose(wave, w)

    parameters = np.concatenate(parameters)
    spectra = np.concatenate(spectra)

    with h5py.File(outname, "w") as out:
        w = out.create_dataset("wavelengths", data=wave)
        p = out.create_dataset("parameters", data=parameters)
        s = out.create_dataset("spectra", data=spectra)