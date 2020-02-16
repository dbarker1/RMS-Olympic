import h5py

files = ["analysis", "run_parameters", "snapshots"]
for file_name in files:
    ana = h5py.File('raw_data/' + file_name + '/' + file_name + '_test1.h5', 'r')
    print("Opened file " + file_name)
    for key in list(ana.keys()):
        print(key)
        for obj in list(ana[key]):
            print(obj)
        #print("The shape of " + key + " is " + str(ana[key].shape) + ".")