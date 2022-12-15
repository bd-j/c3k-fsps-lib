#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Make the low resolution flux files for the C3K library.
"""

import time
from functools import partial
from itertools import product

from multiprocessing import Pool

from ckc.make_full_h5 import full_h5

if __name__ == "__main__":

    from ckc.utils import get_ckc_parser
    # key arguments are:
    #  * --feh
    #  * --afe
    #  * --np
    #  * --ck_vers
    #  * --spec_type
    #  * --fulldir
    #  * --basedir
    parser = get_ckc_parser()
    args = parser.parse_args()
    args.fulldir = args.fulldir.format(args.ck_vers)

    ncpu = args.np
    if ncpu == 1:
        M = map
    else:
        pool = Pool(ncpu)
        M = pool.map

    # --- Metallicities to loop over/map ---
    #fehlist = full_params['feh']
    if args.feh < -10:
        fehlist = [-4.0, -3.5, -3.0, -2.75, -2.5, -2.25, -2.0,
                   -1.75, -1.5, -1.25, -1.0,
                   -0.75, -0.5, -0.25, 0.0, 0.25, 0.5]
    else:
        fehlist = [args.feh]
    if args.afe < -10:
        afelist = [-0.2, 0.0, 0.2, 0.4, 0.6]
    else:
        afelist = [args.afe]

    metlist = list(product(fehlist, afelist))
    print(len(metlist))
    task = partial(full_h5, args=args)

    ts = time.time()
    filenames = list(M(full_h5, list(metlist)))
    dur = time.time() - ts

    print(filenames)
    print('took {}s'.format(dur))
    try:
        pool.terminate()
    except:
        pass
