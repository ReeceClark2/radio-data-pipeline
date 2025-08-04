from astropy.io import fits
import numpy as np
import sys

def compare_headers(header1, header2, ext_name="PRIMARY"):
    keys1 = set(header1.keys())
    keys2 = set(header2.keys())

    all_keys = sorted(keys1.union(keys2))
    print(f"\nComparing headers in extension: {ext_name}")
    for key in all_keys:
        if key not in header1:
            print(f"  - Key '{key}' only in second file.")
        elif key not in header2:
            print(f"  - Key '{key}' only in first file.")
        elif header1[key] != header2[key]:
            print(f"  - Key '{key}' differs:")
            print(f"    First:  {header1[key]}")
            print(f"    Second: {header2[key]}")

def compare_data(data1, data2, ext_name="PRIMARY"):
    if data1 is None and data2 is None:
        return

    print(f"\nComparing data in extension: {ext_name}")

    if data1 is None or data2 is None:
        print("  - One file has data, the other doesn't.")
        return

    if data1.shape != data2.shape:
        print(f"  - Data shapes differ: {data1.shape} vs {data2.shape}")
        return

    # Check if the data is a record array (i.e., table with named columns)
    if hasattr(data1, "dtype") and data1.dtype.names:
        col_names = data1.dtype.names
        for col in col_names:
            col1 = data1[col]
            col2 = data2[col]

            # Choose comparison method
            if np.issubdtype(col1.dtype, np.floating):
                diff_mask = ~np.isclose(col1, col2, equal_nan=True)
            else:
                diff_mask = col1 != col2

            differing_indices = np.where(diff_mask)[0]
            for idx in differing_indices:
                # print(f"  - Row {idx}, Column '{col}': First = {col1[idx]}, Second = {col2[idx]}")
                pass
    else:
        # Fallback for image-like or non-column data
        if np.issubdtype(data1.dtype, np.floating):
            diff = ~np.isclose(data1, data2, equal_nan=True)
        else:
            diff = data1 != data2

        indices = np.argwhere(diff)
        if indices.size == 0:
            print("  - Data is identical.")
            return

        print(f"  - Data differs in {len(indices)} elements:")
        for idx in indices:
            idx_tuple = tuple(idx)
            val1 = data1[idx_tuple]
            val2 = data2[idx_tuple]
            print(f"    At index {idx_tuple}: First = {val1}, Second = {val2}")


def compare_fits(file1, file2):
    with fits.open(file1) as hdul1, fits.open(file2) as hdul2:
        n_ext = max(len(hdul1), len(hdul2))
        print(f"Comparing {n_ext} extensions...")

        for i in range(n_ext):
            if i >= len(hdul1):
                print(f"\nExtension {i}: only in second file.")
                continue
            if i >= len(hdul2):
                print(f"\nExtension {i}: only in first file.")
                continue

            extname = hdul1[i].name if hdul1[i].name else f"EXT{i}"
            compare_headers(hdul1[i].header, hdul2[i].header, extname)
            compare_data(hdul1[i].data, hdul2[i].data, extname)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_fits.py file1.fits file2.fits")
        sys.exit(1)

    file1, file2 = sys.argv[1], sys.argv[2]
    compare_fits(file1, file2)
