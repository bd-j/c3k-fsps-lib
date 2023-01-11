#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from itertools import product

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages
from read_isoc import isoc_table

import h5py

sps_home = os.environ["SPS_HOME"]
template = "{}/{}_feh{:+3.2f}_afe{:+2.1f}.{}.h5"


def isoc_table(feh, afe, itype="MIST"):

    if feh >= 0:
        fsign = "p"
    else:
        fsign = "m"
    if afe >= 0:
        asign = "p"
    else:
        asign = "m"

    if feh == 0.5:
        feh = 0.4
    if itype == "MIST":
        ifile = f"isoc_MIST2_z{sign}{np.abs(feh):3.2f}_afe{afe:+2.1f}.dat"
        absfile = os.path.join(sps_home, "ISOCHRONES", itype, ifile)
    else:
        absfile = f"../data/isoc/mist2_iso/isoc_feh_{fsign}{np.abs(feh*100):03.0f}_afe_{asign}{np.abs(afe*10):01.0f}_vvcrit0.4_full.dat"
    with open(absfile, "r") as f:
        # drop the comment hash and mags field
        header = f.readline().split()[1:]
        #header += list_filters()
    cmd_data = np.loadtxt(absfile, comments="#", dtype=np.dtype([(n, float) for n in header]))
    return cmd_data, absfile


if __name__ == "__main__":

    sedname, prefix = "lr", "c3k_lr"
    sdir = f"../output/{sedname}"
    idir  = os.path.join(sps_home, "ISOCHRONES", "MIST")

    afelist = [-0.2, 0.0, 0.2, 0.4, 0.6]
    fehlist = np.genfromtxt(os.path.join(sdir, "for_fsps", f"{prefix}_zlegend.dat"))
    metlist = list(product(fehlist, afelist))

    pdf = PdfPages("c3k_hrd.pdf")

    for feh, afe in metlist:
        #feh, afe = -0.25, 0.2
        sname = template.format(sdir, "c3k_v1.3", feh, afe, sedname)
        with h5py.File(sname.replace(".h5", ".rect.h5"), "r") as sdat:
            params = sdat["parameters"][:]
            if "extended" in sdat:
                ext = sdat["extended"][:]
            else:
                ext = np.zeros(len(params), dtype=bool)

        fig, ax = pl.subplots(figsize=(9, 6))
        ax.plot(params["logt"], params["logg"], linestyle="", marker="o", color="gray", label="extrapolated")
        ax.plot(params["logt"][ext == 0], params["logg"][ext == 0], linestyle="", marker="o", color="black", label="native")

        try:
            isoc, _ = isoc_table(feh, afe, itype="bleeding edge")
            for logage in [6, 7, 8, 9, 10]:
                sel = np.isclose(isoc["log(age)"], logage, 1e-3)
                tisoc = isoc[sel]
                tisoc = tisoc[np.argsort(tisoc["Mini"])]
                ax.plot(tisoc["logt"], tisoc["logg"], label=f"MIST2: log(Age)={logage:.2f}")
        except(FileNotFoundError):
            pass

        ax.set_xlim(4.79, 3.31)
        ax.set_ylim(5.8, -1.2)
        ax.set_xlabel(r"$\log {\rm T}_{eff}$", fontsize=18)
        ax.set_ylabel(r"$\log g$", fontsize=18)
        ax.legend()
        ti = rf"$[{{\rm Fe}}/{{\rm H}}]={feh:+3.2f}, [\alpha/{{\rm Fe}}]={afe:+2.1f}$"
        ax.set_title(ti, fontsize=18)
        pdf.savefig(fig)
        pl.close(fig)
    pdf.close()