import numpy as np
from astropy.io import fits

with fits.open("../data/cyg_a_low_hI/0137048.fits") as orig, fits.open("../data/cyg_a_low_hI_mod/0137048.fits") as mod:
    print("Unique CRVAL1 values (original):", np.unique(orig[1].data['CRVAL1']))
    print("Unique CRVAL1 values (modified):", np.unique(mod[1].data['CRVAL1']))
    
    print("Unique OBSFREQ values (original):", np.unique(orig[1].data['OBSFREQ']))
    print("Unique OBSFREQ values (modified):", np.unique(mod[1].data['OBSFREQ']))
