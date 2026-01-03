import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.time import Time
import utils
from typing import Sequence


class Merge:
    def __init__(self, file_paths: Sequence[str]):
        '''
        Initialize all files to be merged.
        '''

        if len(file_paths) < 2:
            raise ValueError("At least two FITS files are required.")

        self.file_paths = list(file_paths)
        self.headers = []
        self.tables = []

        for path in self.file_paths:
            with fits.open(path) as hdul:
                hdul.verify("exception")
                self.headers.append(hdul[0].header)
                self.tables.append(Table(hdul[1].data))

        self._validate_tables()

    def _validate_tables(self):
        '''
        Validate all data tables.
        '''

        ref = self.tables[0]

        for i, tbl in enumerate(self.tables[1:], start=1):
            if tbl.colnames != ref.colnames:
                raise ValueError(f"Column mismatch in file {i}.")

            for name in ref.colnames:
                if tbl[name].dtype != ref[name].dtype:
                    raise ValueError(f"Dtype mismatch in column '{name}' (file {i}).")

    def merge(self):
        '''
        Merge all data sections using vstack and maintain 
        the first file header as the primary header.
        '''

        merged = vstack(self.tables, metadata_conflicts="silent")
        header = self.headers[0]

        utils.save(self.file_paths[0], header, merged, "merge")


if __name__ == "__main__":
    filepaths = [
        "C:/Users/starb/Downloads/Raw/0144717daisy.fits",
        "C:/Users/starb/Downloads/Raw/0144718map.fits",
    ]

    m = Merge(filepaths)
    m.merge()
