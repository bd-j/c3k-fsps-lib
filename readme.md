# c3k-lib

Code here is for generating low-resolution spectra files for (python-)fsps,
using the c3k library.

## Resolution

Resolution is specified in a yml file as tuples of `(lambda_low, lambda_hi, R, usefft)`

Note the c3k_v1.3 default resolutions are

* 0.009micron - 0.0912micron; opacity binning with wavelength spacing lambda/Delta_lambda=4500
* 0.09121micron - 2.5micron; opacity binning with wavelength spacing lambda/Delta_lambda=100000
* 2.5micron - 40micron; opacity binning with wavelength spacing lambda/Delta_lambda=4500

Beyond 40 microns we stitch a Rayleigh-Jeans tail onto the spectra.

## Procedure

1. Install.  You should have FSPS already installed and `$SPS_HOME` environment
   variable set.

   ```sh
   cd $SCRATCH/conroy_lab/$USER
   git clone git@github.com:bd-j/c3k-fsps-lib.git
   cd c3k-fsps-lib
   module purge
   module load Anaconda3/2020.11
   conda env create -f environment.yml
   source activate c3k
   cd ..
   ```

2. Check that all the relevant flux files are there

   ```sh
   fdir=/n/holystore01/LABS/conroy_lab/Lab/bdjohnson/data/kurucz/c3k_v1.3/fullres
   for f in $fdir/*full*h5; do ls ${f/.full./.flux.}; done
   ```

   If they are not, then you need to make them using `make_flux.py`, driven by
   `ody_flux.sh`

3. Make the low resolution SED H5 files, and the H5 files of the SEDs
   interpolated to the FSPS logt-logg gridpoints (based on BaSeL or, in new
   versions, c3k).  Also make binary format versions for FSPS itself, and
   metadata files (wavelength, resolution, zlegend)

   ```sh
   cd jobs/
   libname=nirspec # name of the segments_<libname>.yml file
   sbatch --export=ALL,libname=${libname} --array=0-64 ody_resample.sh
   ```

   Note that you may wish to change the resampling here, the metallicities to
   consider, or the output filename and directory.  The output names are given
   by the supplied `libname` and the metallicities to consider are
   specified in `c3k_resample.py`.  The number of jobs in the array should be
   equal to the number of feh-afe pairs (usually 13 * N_afe)

   There are now several sets of binary files and ancillary files in the
   `output/${libname}/for_fsps/` directory (by default) that can be moved to the
   ```$SPS_HOME/SPECTRA/C3K``` directory.  By altering `sps_vars.f90` you can
   choose to use these spectra.

4. If you're me the `output/${libname}` directory should be copied to
   `/n/holystore01/LABS/conroy_lab/Lab/bdjohnson/data/kurucz/${ck_vers}/fsps-lib/`

   ```sh
   cd jobs/
   libname=nirspec # name of the segments_<libname>.yml file
   sbatch --export=ALL,libname=${libname} ody_copy.sh
   ```

5. Implement in FSPS

   This depends a bit on what version of FSPS you have. Instructions her are for
   the 'older' FSPS. Assuming you have downloaded FSPS, the first thing to do is
   copy the relevant files to the FSPS repo.  Note this changes version tracked
   files!  The prefix also changes.

   ```sh
   cp <path/to/output/libname>/for_fsps/<prefix>*bin $SPS_HOME/SPECTRA/C3K/
   cp <path/to/output/libname>/for_fsps/<prefix>*dat $SPS_HOME/SPECTRA/C3K/
   cp <path/to/output/libname>/for_fsps/<prefix>_zlegend.dat $SPS_HOME/SPECTRA/C3K/zlegend.dat
   cp <path/to/output/libname>/for_fsps/<prefix>.wave $SPS_HOME/SPECTRA/C3K/<prefix>.lambda
   ```

   Then, you need to change the `sps_vars.f90` code.  You'll want to specify the
   the C3K library using the pre-compiler directives, and you'll need to change
   environment variables in the ``elif (C3K)`` block; the values for `nzinit`
   and `nspec` can be obtained by `wc` on the `*.lambda` and `*_zlegend.dat`
   files, and the `spec_type` variable (and its length!) could be changed to
   e.g. `c3k_ns`. Here's what the diff of `sps_vars.f90` looks like with
   `nz=11`, `nspec=13749` and a prefix of `c3k_ns`:

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
   +        CHARACTER(6), PARAMETER :: spec_type = 'c3k_ns'
   +        INTEGER, PARAMETER :: nzinit=11
   +        INTEGER, PARAMETER :: nspec=13749
   ```

   You can now recompile fsps (`cd $SPS_HOME/src; make clean; make all`) and it
   should use the new C3K spectral library.

6. Implement in python-fsps

   This requires a [development install](https://dfm.io/python-fsps/current/installation/#installing-development-version)
   of python-fsps. You need to copy the files to `$SPS_HOME` as above. Then, after
   cloning the repo, you'll need to edit `fsps/src/fsps/libfsps/src/sps_vars.f90`
   in the same way as described above. Finally, install with the C3K library selected.
   This looks like:

   ```sh
   cp <path/to/output/libname>/for_fsps/c3k_ns* $SPS_HOME/SPECTRA/C3K/
   cp <path/to/output/libname>/for_fsps/c3k_ns_zlegend.dat $SPS_HOME/SPECTRA/C3K/zlegend.dat

   git clone --recursive https://github.com/dfm/python-fsps.git
   <change python-fsps/fsps/src/fsps/libfsps/src/sps_vars.f90 as in step 5>

   cd python-fsps
   python -m pip uninstall fsps
   FFLAGS="-DMILES=0 -DC3K=1" python -m pip install .
   ```

7. Test python-fsps implementation

   You can test that this worked by running a script in the test/ directory and
   looking at the output.  (Make sure to do this the environment where you
   installed python-fsps above) . This may take some time, as the SSPs are
   regenerated several times.

   ```sh
   cd c3k-fsps/lib/tests
   python fsps_feature_demo.py
   open features.pdf
   ```

   You can also plot the grid for every metallicity and the corresponding
   isochrones with `make_hrd.py`

   It might also be instructive to compare the SSP spectra to another library (e.g. MILES)