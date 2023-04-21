import argparse

import h5py

from libtmdpdf.utils.print_h5 import print_h5

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Print the structure of HDF5 files.")
    parser.add_argument(
        "files", metavar="F", type=str, nargs="+", help="an HDF5 file path"
    )
    args = parser.parse_args()

    for file in args.files:
        try:
            print(f"{file}:")
            with h5py.File(file, "r") as f:
                print_h5(f, 1)
        except:
            print("\tERROR")
