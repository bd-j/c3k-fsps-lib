from itertools import product
import os, glob
import numpy as np
import h5py

cc_name = "/n/holystore01/LABS/conroy_lab/Lab/cconroy/kurucz/grids/c3k_v2.3"
ext = "spec"
dirname = "/n/holystore01/LABS/conroy_lab/Lab/bdjohnson/data/kurucz/c3k_v2.3/optfal/spec"

fehlist = [-2.0, -1.75, -1.5, -1.25, -1.0,
           -0.75, -0.5, -0.25, 0.0, 0.25, 0.5,
           -2.25, -2.5, -2.75, -3.0, -3.5, -4.0]
afelist = [-0.2, +0.0, +0.2, +0.4, +0.6]
metlist = list(product(fehlist, afelist))

for feh, afe in metlist:
    zstring = f"feh{feh:+3.2f}_afe{afe:+2.1f}"
    libname = f"{dirname}/c3k_v2.3_{zstring}.spec.h5"
    c3k = h5py.File(libname, "r+")
    params = c3k["parameters"][:]
    c3k_wave = c3k["wavelengths"][:]
    print(f"operating on {zstring}")

    for i, p in enumerate(params):
        t = 10**p["logt"]
        g = f"{p['logg']:4.2f}"

        kstring = f"t{t:05.0f}g{g:.4s}"
        fn = f"{cc_name}/at12_{zstring}/spec_optfal/at12_{zstring}_{kstring}.{ext}"
        if not os.path.exists(fn):
            print(f"could not find {fn}")
            continue

        fulldf = np.loadtxt(fn)
        nc = fulldf.shape[-1]
        assert nc == 3

        wave = np.array(fulldf[:, 0])
        full_spec = np.array(fulldf[:, 1])  # spectra
        full_cont = np.array(fulldf[:, 2])  # continuum

        new = np.interp(c3k_wave, wave, full_spec, left=0.0, right=0.0)
        sel = new > 0
        assert sel.sum() > 0, f"no spectra at {i}: {fn}"
        c3k["spectra"][i, sel] = new[sel]

        new = np.interp(c3k_wave, wave, full_cont, left=0.0, right=0.0)
        sel = new > 0
        assert sel.sum() > 0, f"no continuum at {i}: {fn}"
        c3k["continuua"][i, sel] = new[sel]

    c3k.close()
