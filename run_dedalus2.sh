#!/usr/bin/env bash

name=$1
cores=$2

git clone https://github.com/dbarker1/rm-history.txt
mv rm-history.txt $name
cd $name
rm notes_RBconvection_script.pdf
cp ~/run_param_file2.py .
git checkout origin/Corrected_Rotation

mpiexec -np $cores python3 $(pwd)/anelastic_RB.py
rm -r __pycache__
cp -r raw_data raw_data_cp

merge.py raw_data/snapshots --cleanup
merge.py raw_data/analysis --cleanup
merge.py raw_data/run_parameters --cleanup

merge_single.py $name
plotting_snapshots.py $name
