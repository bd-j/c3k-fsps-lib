import numpy as np

__all__ = ["ckms", "sigma_to_fwhm", "tiny_number",
           "param_order"]

# Useful constants
ckms = 2.998e5
sigma_to_fwhm = 2 * np.sqrt(2 * np.log(2))
tiny_number = 1e-33


param_order = ['t', 'g', 'feh', 'afe']