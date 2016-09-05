#!/bin/sh


dir_out=/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20150111_LB_Sensitivity/
aperture_diameter_mm=400
T_bath=0.1
python Mapping_speed_NdivNET2.py ${dir_out}/data $aperture_diameter_mm $T_bath

dir_summary=$dir_out/out
python mapping_speed_summary_table.py ${dir_out}/data $dir_summary ${aperture_diameter_mm} 
