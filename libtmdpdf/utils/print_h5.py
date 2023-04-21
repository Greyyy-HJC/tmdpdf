import argparse

import h5py


def print_h5(group, level=0):
    if not isinstance(group, h5py.Group):
        print(level * "\t", group)
        return
    else:
        for key in group.keys():
            print(level * "\t" + key + ":")
            subgroup = group[key]
            print_h5(subgroup, level + 1)


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
