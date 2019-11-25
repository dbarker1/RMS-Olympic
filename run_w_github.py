#!/usr/bin/env python3

import subprocess as sp
import sys
path_to_sim = "/home/djb236/dedalus_sims"
import os

if (len(sys.argv) < 2):
	print("Must provide number of cores")
	exit(1)

CORES = sys.argv[1]

class Param:
    def __init__(self, name_in, val_in_dot):
        self.name = name_in[:-1]
        self.val = val_in_dot
        if (val_in_dot.find(".") != -1):
            self.val_p = change_to_ps(self.val)

    def change_param(self, val_in):
        self.val = val_in
        if (val_in.find(".") != -1):
            self.val_p = change_to_ps(self.val)

params = []

def change_to_ps(string):
    name = ""
    for char in string:
        if (char == "."):
            name += "p"
        else:
            name += char
    return name

def read_params(string):
    for line in string.split("\n"):
        if (line != ""):
            if (line[0] != '#' and line[:6] != "import" and line != ""):
                if (line.split("=")[0].find("theta") != -1):
                    params.append(Param("theta ", "1 - np.exp(-Np/m)"))
                else:
                    params.append(Param(line.split("=")[0], line[line.find("="):].split(" ")[1]))

def get_val(param_str, print_dot):
    for param in params:
        if (param.name == param_str):
            if (print_dot):
                return param.val
            else:
                if (param.val.find("np.pi") != -1):
                      return "pi_" + param[len("np.pi")+1:]
                return change_to_ps(param.val)

def edit_param(name):
    for param in params:
        if (name == param.name):
            param.change_param(input("New value for " + param.name + ": "))

def print_all():
    out = ""
    for param in params:
        out += param.name + " = " + param.val + "\n"
    return out

def print_msg(msg):
    stars = ""
    for char in range(len(msg) + 2):
        stars += "*"
    print("\n\n" + stars)
    print(" " + msg + " ")
    print(stars)

f_op = open("/home/djb236/RMS-Olympic/run_param_file2.py", "r")
read_params(f_op.read())
f_op.close()

print("Initial parameters:")
print(print_all())
edit_param(input("Parameter to change: "))
print("\nNew parameters:")
print(print_all())
param_path = "Np_" + get_val("Np", False) + "/Ra_" + get_val("Ra", False) + "/Ta_" + get_val("Ta", False) + "Lat_" + get_val("phi", False)

rot_no = 0
for item in os.listdir(path_to_sim):
    if (item.find("rot") != -1):
        if (float(item.split("rot")[1]) > rot_no):
            rot_no = int(item.split("rot")[1])

name = "rot" + str(rot_no+1)
sp.call("mkdir " + name, shell=True)
print_msg("Getting anelastic script from Github...")
sp.call("cd " + name + " && wget https://raw.githubusercontent.com/dbarker1/RMS-Olympic/2.5D_Rotation/anelastic_RB.py", shell=True)

cwd = path_to_sim + "/" + name
f_wr_op = open(cwd + "/run_param_file2.py", "w")
f_wr_op.write("import numpy as np\n\n" + print_all())
f_wr_op.close()

print_msg("Running dedalus...")
sp.call("cd " + cwd + " && mpiexec -np " + str(CORES) + " python3 " + cwd + "/anelastic_RB.py", shell=True)
print_msg("Removing cache folder...")
sp.call("cd " + cwd + " && rm -r __pycache__ && cp -r raw_data raw_data_cp", shell=True)

print_msg("merge.py:")
sp.call("merge.py " + cwd + "/raw_data/snapshots --cleanup", shell=True)
sp.call("merge.py " + cwd + "/raw_data/analysis --cleanup", shell=True)
sp.call("merge.py " + cwd + "/raw_data/run_parameters --cleanup", shell=True)

print_msg("merge_single.py:")
sp.call("cd " + cwd + " && merge_single.py " + name, shell=True)
print_msg("plotting_snapshots.py:")
sp.call("cd " + cwd + " && plotting_snapshots.py " + name, shell=True)

sp.call("mkdir -p ~/RMS-Olympic/RESULTS/" + param_path, shell=True)
sp.call("cd " + cwd + "/" + name + "_figs && cp -r * ~/RMS-Olympic/RESULTS/" + param_path, shell=True)
print_msg("Pushing to Github...")
sp.call("cd ~/RMS-Olympic && git add . && git commit -m 'NEW: " + name + "' && git push", shell=True)
print_msg("Done")
