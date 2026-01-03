import os
from astropy.io import fits
from astropy.table import Table

class Corrupt:
    def __init__(self, file_path: str):
        self.filepath = file_path
        self.missing_values = []

        with fits.open(self.filepath) as hdul:
            self.header = hdul[0].header.copy()
            print(repr(self.header))
            self.data = Table(hdul[1].data)

    def corrupt(self):
        if 'END' in self.header:
            del self.header['END']
            self.missing_values.append('END')

    def save(self, output_path=None):
        if output_path is None:
            base, ext = os.path.splitext(self.filepath)
            output_path = f"{base}_corrupted{ext}"

        primary_hdu = fits.PrimaryHDU(header=self.header)

        if 'END' in primary_hdu.header and 'END' in self.missing_values:
            del primary_hdu.header['END']
            # pass

        table_hdu = fits.BinTableHDU(data=self.data)
        hdulist = fits.HDUList([primary_hdu, table_hdu])

        hdulist.writeto(output_path, overwrite=True, output_verify='ignore')

if __name__ == "__main__":
    filepath = "C:/Users/starb/Downloads/0136873(1).fits"
    c = Corrupt(filepath)
    c.corrupt()
    c.save()