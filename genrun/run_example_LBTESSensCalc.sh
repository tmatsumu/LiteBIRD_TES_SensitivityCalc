#!/bin/sh

# Change the following path, dir_out, to specify your output directory
dir=/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20160906_LiteBIRDTESSensitivityCalc/

python $dir/src/Mapping_speed_NdivNET2.py ${dir}/data/

# order of the code execusion is 
# given global_par.py
# 1) Mapping_speed_NdivNET2.py
# exe_mappingspeed.py
# lib_mappingspeed.py
