from astropy.io import fits

# Open the FITS file
with fits.open("C:/Users/starb/Downloads/0136993.fits", mode="update") as hdul:
    header = hdul[0].header  # Primary header

    # Modify BMAJ and BMIN
    header['BMAJ'] = 0.66537725863748  # New value for BMAJ
    header['BMIN'] = 0.66537725863748  # New value for BMIN

    # Optionally save as new file
    hdul.writeto("output_modified.fits", overwrite=True)

print("FITS file modified and saved as '0136993_output_modified.fits'")
