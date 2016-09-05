#!/bin/sh

dir_in=/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20151223_LB_HFT_Sensitivity/

dir_src=${dir_in}/src/
#dir_out=${dir_in}//data/HWP0p01diffTrans/
dir_out=${dir_in}//data/LoadingVariation_DipoleAmplitude/
mkdir -p ${dir_out}
python ${dir_src}/main_CalcLoading_v02.py ${dir_out}

# order of the code execusion is 
# given global_par.py
# 1) Mapping_speed_NdivNET2.py
# exe_mappingspeed.py
# lib_mappingspeed.py

#dir_summary=$dir_out/out
#python mapping_speed_summary_table_20141207_Tmir.py ${dir_out}/data $dir_summary
