#! /bin/bash

source /opt/Miniconda/miniconda37/etc/profile.d/conda.sh
conda activate dedalus

#cp run_param_file.py run_param_file_$1.py

mpiexec -n 4 python3 main.py
exit
