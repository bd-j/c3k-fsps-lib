# Output

Place to put separate directories containing hdf5 and binary files for custom
spectral libraries.

## Description

These spectra were created from the native resolution C3K SYNTHE output .spec
files (lambda/Delta_lambda \sim 100,000, \lambda\lambda=900-25,0000 angstrom).
The native SYNTHE resolution files were convolved (via FFT in the middle, and
directly at the edges) with a Gaussian in velocity (or ln(lam)) space to yield a
final output resolution of \lambda/(\Delta\lambda) = R, where \Delta\lambda is
the FWHM of the Gaussian. In other words, the output velocity resolution is
c/R/2.355 = (10000/R) * 12.7 km/s (dispersion, not FWHM).

Shortward of 1000AA and longward of 2.5micron the spectra are supplemented with
lower-resolution SEDs based on the SYNTHE "flux" output, which is constructed
from total estimated line opacity in bins of wavelength  (with
lambda/Delta_lambda ~ 4500), but it is not based on true spectral synthesis.
The resolution is not well defined, and these "flux" outputs are further
smoothed.

## Files

Located in `output/<libname>` where `<libname>` is the name of the library.  The
_segments_ file that specifies the resolutions is also at this location.

### .h5 Files

These are the c3k grid at the resolution indicatedThe HDF5 data are arranged in three HDF5 datasets:

1. `wavelengths`: A vector giving the restframe vacuum wavelength scale, in angstroms.
   It has N_wave elements.

2. `spectra`: A 2 dimensional array of f\nu spectra, _not_ continuum
   normalized, in units of erg/s/cm^2/Hz/sr. Multiply by (4\pi R_{star})^2 to obtain
   absolute units (erg/s/Hz), where L_{star} = 4\pi\sigma_{SB} R_{star}^2 T_{star}^4. The array
   has a shape of N_models x N_wave.

3. `parameters`: A structured array with N_stars rows giving the stellar
   parameters of each spectrum. The 4 fields or columns of the structured array
   are:
      * `logt` (log_10 of the effective temperature, K),
      * `logg` (log_10 of the gravity),
      * `feh` (the value of [Fe/H], i.e. log_10(n_{Fe}/n_H) - log_10(n_{Fe,\odot}/n_{H, \odot})) and
      * `afe` the value of [\alpha/Fe] (always 0.0)
   Note that for C3K the solar metallicity is Z=0.0134, and the solar abundances
   are from Asplund 2009.

### fsps.h5 Files

Similar to above, but the spectra have been interpolated onto the grid of
logTeff and logg appropriate for BaSeL and FSPS.  The `parameters` dataset only
includes `logt` and `logg`.

### for_fsps directory

This directory contains the files appropriate for FSPS.  All the binary spectral
files are prefixed `c3k_alpha+0.0`.  There is also a set of *zlegend.dat and
*.lambda files and resolution file.

## Implement in FSPS

Assuming you have downloaded FSPS, the first thing to do is copy the relevant
files to the FSPS repo.  Note this changes version tracked files!  In principle
you could change the prefix, but that's more complicated.
```sh
cp <path/to/output/libname>/for_fsps/c3k_afe+0.0* $SPS_HOME/SPECTRA/C3K/
cp <path/to/output/libname>/for_fsps/c3k_afe+0.0_zlegend.dat $SPS_HOME/SPECTRA/C3K/zlegend.dat
```

You'll need to change environment variables in the ``elif (C3K)`` block; the
values for `nzinit` and `nspec` can be obtained by `wc` on the `*.lambda` and
`*_zlegend.dat` files, and the `spec_type` variable (and its length!) could be
changed to e.g. `c3k_ns_afe+0.0`. Finally, the appropriate zlegend file should
be copied to `zlegend.dat` (removing the prefix).  Here's what the diff of
`sps_vars.f90` looks like with nz=11 and nspec=13749 and no change to the prefix:

```diff
@@ 7
-#define MILES 1
+#define MILES 0

@@ 16
-#define C3K 0
+#define C3K 1

   @@ 239
   #elif (C3K)
      REAL(SP), PARAMETER :: zsol_spec = 0.0134
-        CHARACTER(11), PARAMETER :: spec_type = 'c3k_afe+0.0'
-        INTEGER, PARAMETER :: nzinit=11
-        INTEGER, PARAMETER :: nspec=11149  !46666 !47378 !, 26500
+        CHARACTER(11), PARAMETER :: spec_type = 'c3k_afe+0.0'
+        INTEGER, PARAMETER :: nzinit=11
+        INTEGER, PARAMETER :: nspec=13749  !46666 !47378 !, 26500
   ```

You can now recompile fsps (`cd $SPS_HOME/src; make clean; make all`) and it
should use the new C3K spectral library

## Implement in python-fsps

This requires a [development install]
(https://dfm.io/python-fsps/current/installation/#installing-development-version)
of python-fsps. After cloning the repo, you'll need to edit
`fsps/src/fsps/libfsps/src/sps_vars.f90` in the same way as described above.
Then, install with the C3K library selected.  This looks like:

```sh
git clone --recursive https://github.com/dfm/python-fsps.git
<edit python-fsps/fsps/src/fsps/libfsps/src/sps_vars.f90>
cd python-fsps
pythn -m pip uninstall fsps
FFLAGS="-DMILES=0 -DC3K=1" python -m pip install .
```