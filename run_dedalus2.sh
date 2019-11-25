#!/usr/bin/env bash

name=$1
cores=$2

mkdir $name
cd $name
wget https://raw.githubusercontent.com/dbarker1/rm-history.txt/2.5D_Rotation/anelastic_RB.py
cp ~/run_param_file2.py .

mpiexec -np $cores python3 $(pwd)/anelastic_RB.py
rm -r __pycache__
cp -r raw_data raw_data_cp

merge.py raw_data/snapshots --cleanup
merge.py raw_data/analysis --cleanup
merge.py raw_data/run_parameters --cleanup

merge_single.py $name
plotting_snapshots.py $name
