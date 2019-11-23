#!/usr/bin/env python3

import subprocess as sp
path_to_sim = "/home/djb236/dedalus_sims"
CORES = 10
import os

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
                params.append(Param(line.split("=")[0], line[line.find("="):].split(" ")[1]))

def get_val(param_str, print_dot):
    for param in params:
        if (param.name == param_str):
            if (print_dot):
                return param.val
            else:
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

f_op = open("/home/djb236/run_param_file2.py", "r")
f_str = f_op.read()
f_op.close()

read_params(f_str)
print("Initial parameters:")
print(print_all())
edit_param(input("Parameter to change: "))
print("\nNew parameters:")
print(print_all())
param_path = "Np_" + get_val("Np",False) + "/Ra_" + get_val("Ra", False) + "/Ta_" + get_val("Ta", False)
print(param_path)
rot_no = 0
for item in os.listdir(path_to_sim):
    if (item.find("rot") != -1):
        if (float(item.split("rot")[1]) > rot_no):
            rot_no = int(item.split("rot")[1])

name = "rot" + str(rot_no+1)
sp.call("mkdir " + name, shell=True)
sp.call("cd " + name + " && wget https://raw.githubusercontent.com/dbarker1/rm-history.txt/2.5D_Rotation/anelastic_RB.py", shell=True)

cwd = path_to_sim + "/" + name
f_wr_op = open(cwd + "/run_param_file2.py", "w")
f_wr_op.write("import numpy as np\n\n" + print_all())
f_wr_op.close()

sp.call("mpiexec -np " + str(CORES) + " python3 " + cwd + "/anelastic_RB.py", shell=True)
sp.call("cd " + cwd + " && rm -r __pycache__ && cp -r raw_data raw_data_cp", shell=True)

sp.call("merge.py " + cwd + "/raw_data/snapshots --cleanup", shell=True)
sp.call("merge.py " + cwd + "/raw_data/analysis --cleanup", shell=True)
sp.call("merge.py " + cwd + "/raw_data/run_parameters --cleanup", shell=True)

sp.call("merge_single.py " + name, shell=True)
sp.call("plotting_snapshots.py " + name, shell=True)

sp.call("mkdir -p ~/rm-history.txt/RESULTS" + param_path, shell=True)
sp.call("cd " + cwd + "/" + name + "_figs && cp -r * ~/rm-history.txt/RESULTS/", shell=True)
sp.call("cd ~/rm-history.txt && git add . && git commit -m '" + name + "' && git push", shell=True)
