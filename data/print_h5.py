#!python
# %%
import h5py
import pandas as pd
import sys
# %%
# %%
def print_h5(group,level=0):
    if(isinstance(group,h5py.Group)==False):
        print(level*'\t',group)
        return
    else:
        for key in group.keys():
            print(level*'\t'+key+':')
            subgroup=group[key]
            print_h5(subgroup,level+1)

# %%
if __name__=='__main__':
    if len(sys.argv)==1 or \
        (len(sys.argv)==2 and (sys.argv[1]=='-h' or sys.argv[1]=='--help')):
        print("""usage: \
to print the structure of HDF5 files.
\tshell>> python print_h5.py file1 [file2 ...]""")
        exit(0)

    for file in sys.argv[1:]:
        try:
            print(file+':')
            print_h5(h5py.File(file,'r'),1)
            pass
        except Exception:
            print('\tERROR')
            pass
        