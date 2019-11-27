#!/usr/bin/env python3
"""
Perform after execution of the merge.py script to combine all .h5 data files
into a singular .h5 file.

WARNING: This script may produce an extremely large output file.
"""

import pathlib
from dedalus.tools import post
import sys

# direc is the folder in which the various data folders are located

#-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#
############################# These strings need editting #############################
#-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#
#run_name = "test1"
run_name = sys.argv[1]
direc_folder= "raw_data/"
#-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#
#######################################################################################
#-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#   #-----#


subfolder = "snapshots"
direc = direc_folder + subfolder
print(direc_folder)
# name of the resulting output file
output_name = subfolder + "_" + run_name

set_paths = list(pathlib.Path(direc).glob(subfolder + "_s*.h5"))
post.merge_sets(direc + "/" + output_name + ".h5", set_paths, cleanup=True)

# direc is the folder in which the various data folders are located
subfolder = "analysis"
direc = direc_folder + subfolder
# name of the resulting output file
output_name = subfolder + "_" + run_name

set_paths = list(pathlib.Path(direc).glob(subfolder + "_s*.h5"))
post.merge_sets(direc + "/" + output_name + ".h5", set_paths, cleanup=True)

# direc is the folder in which the various data folders are located
subfolder = "run_parameters"
direc = direc_folder + subfolder
# name of the resulting output file
output_name = subfolder + "_" + run_name

set_paths = list(pathlib.Path(direc).glob(subfolder + "_s*.h5"))
post.merge_sets(direc + "/" + output_name + ".h5", set_paths, cleanup=True)
