import h5py
import sys

run_name = sys.argv[1]

files = ["analysis", "run_parameters", "snapshots"]
for file_name in files:
    ana = h5py.File('raw_data/' + file_name + '/' + file_name + '_' + run_name + '.h5', 'r')
    print("\nFile " + file_name + ":")
    for key in list(ana.keys()):
        print(key)
        for obj in list(ana[key]):
            print(ana[key][obj])