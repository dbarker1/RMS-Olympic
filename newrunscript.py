import subprocess as sp

Ra = "1e6"
Tas = ["2e5", "3e5", "6e5", "7e5", "8e5", "9e5"]

store = "/home/djb236/BOUSINESSQ/"

for Ta in Tas:
	dir = store + Ra + "/" + Ta
	sp.call("mkdir -p " + dir, shell=True)
	sp.call("cd " + dir + " && mpiexec -np 10 python3 ~/RMS-Olympic/rot_pseudo2d_anaRB.py 0 " + Ra + " " + Ta + " 30", shell=True)