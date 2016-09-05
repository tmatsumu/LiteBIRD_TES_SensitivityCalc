import numpy as np
import pylab as py
import os 
import sys
import lib_mappingspeed as lib_ms

#######################################################################################################
wafer_num = [8, 8, 8, 5, 5, 5]
#;wafer_num = [1, 1, 1, 1, 1, 1]

#;if fix(total(wafer_num)) eq 6 then subdir_tag = '1wafer/'
#;if fix(total(wafer_num)) eq 39 then subdir_tag = 'HFP5_LFP8_wafers/'

aperture_diameter_mm = 400. #; mm
freq_str=[60,78,100,140,195,280]

num_freq = len(freq_str)

T_bath = 0.1
T_mir_arr = [2., 2.5, 3.0, 3.5, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, \
                 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.]
num_Tmir = len(T_mir_arr)

NEP = np.zeros((num_freq, num_Tmir))
NET = np.zeros((num_freq, num_Tmir))
Pload = np.zeros((num_freq, num_Tmir))
NEP_ph_w = np.zeros((num_freq, num_Tmir))
NEP_th = np.zeros((num_freq, num_Tmir))
NEP_readout = np.zeros((num_freq, num_Tmir))

dir_in = sys.argv[1]
dir_out = sys.argv[2]
#dir_out = '/Users/tomotakematsumura/work/develop/LiteBIRD/NET/CalNET/results/20140717/summary/'
#openw,lun,dir_out+'NEP_Tmir_'+strtrim(fix(aperture_diameter_mm),2)+'mm_4K.txt',/get_lun

for i in range(0,num_Tmir):
    dir = dir_in+'/'+str(int(T_bath*1e3))+'mK_Dapt'+str(int(aperture_diameter_mm))+'mm'
#    '/Users/tomotakematsumura/work/develop/LiteBIRD/NET/CalNET/results/20140717/'+str(int(T_bath*1e3)) \
#        +'mK_Dapt'+str(int(aperture_diameter_mm))+'mm'

    dir_data = dir+'/data'
    dir_out = dir+'/out'
    os.system('mkdir -p '+dir_out)
   
    Tcmb = 2.73
    c = 3.e8
#; beam_ref = 2.481e-8 ; str
#; eff = 0.64
#;   restore, dir_data+'/eff.save'
#;   restore, dir_data+'/emiss.save'
#######################################################################################################
    for j in range(0,num_freq):
        file_in = dir_data+'/mappingspeed_'+str(freq_str[j])+'GHz_F_T'+str(T_mir_arr[i])+'K.npy'
        output = np.load(file_in)
        if freq_str[j] < 110: ind = np.where((output[10] > 0.01799) & (output[10] < 0.01801))
        if freq_str[j] > 110: ind = np.where((output[10] > 0.01199) & (output[10] < 0.01201))
        NEP_ph_w[j,i] = output[19,ind[0]]
        NEP_th[j,i] = output[20,ind[0]]
        NEP_readout[j,i] = output[37,ind[0]]
        NEP[j,i] = output[26,ind[0]]
        NET[j,i] = output[32,ind[0]]
        Pload[j,i] = output[17,ind[0]]

py.figure()
for j in range(0,num_freq):
    if freq_str[j] == 60:
        py.plot(T_mir_arr, Pload[j,:]*1e12,'ob')
        py.plot(T_mir_arr, Pload[j,:]*1e12,'-b', label='60GHz')
    if freq_str[j] == 78:
        py.plot(T_mir_arr, Pload[j,:]*1e12,'og')
        py.plot(T_mir_arr, Pload[j,:]*1e-12,'-g', label='78GHz')
    if freq_str[j] == 100:
        py.plot(T_mir_arr, Pload[j,:]*1e12,'or')
        py.plot(T_mir_arr, Pload[j,:]*1e12,'-r', label='100GHz')
    if freq_str[j] == 140:
        py.plot(T_mir_arr, Pload[j,:]*1e12,'oc')
        py.plot(T_mir_arr, Pload[j,:]*1e12,'-c', label='140GHz')
    if freq_str[j] == 195:
        py.plot(T_mir_arr, Pload[j,:]*1e12,'om')
        py.plot(T_mir_arr, Pload[j,:]*1e12,'-m', label='195GHz')
    if freq_str[j] == 280:
        py.plot(T_mir_arr, Pload[j,:]*1e12,'oy')
        py.plot(T_mir_arr, Pload[j,:]*1e12,'-y', label='280GHz')
py.xlabel('Temperature of optical system [K]')
py.ylabel('Loading [pW]')
py.legend(loc='best')

py.figure()
for j in range(0,num_freq):
    if freq_str[j] == 60:
        py.plot(T_mir_arr, NEP[j,:]*1e18,'ob')
        py.plot(T_mir_arr, NEP[j,:]*1e18,'-b', label='60GHz')
        py.plot(T_mir_arr, NEP_ph_w[j,:]*1e18,'-.b')
        py.plot(T_mir_arr, NEP_th[j,:]*1e18,'--b')
        py.plot(T_mir_arr, NEP_readout[j,:]*1e18,'.b')

    if freq_str[j] == 78:
        py.plot(T_mir_arr, NEP[j,:]*1e18,'og')
        py.plot(T_mir_arr, NEP[j,:]*1e-18,'-g', label='78GHz')
        py.plot(T_mir_arr, NEP_ph_w[j,:]*1e18,'-.b')
        py.plot(T_mir_arr, NEP_th[j,:]*1e18,'--b')
        py.plot(T_mir_arr, NEP_readout[j,:]*1e18,'.b')

    if freq_str[j] == 100:
        py.plot(T_mir_arr, NEP[j,:]*1e18,'or')
        py.plot(T_mir_arr, NEP[j,:]*1e18,'-r', label='100GHz')
        py.plot(T_mir_arr, NEP_ph_w[j,:]*1e18,'-.b')
        py.plot(T_mir_arr, NEP_th[j,:]*1e18,'--b')
        py.plot(T_mir_arr, NEP_readout[j,:]*1e18,'.b')

    if freq_str[j] == 140:
        py.plot(T_mir_arr, NEP[j,:]*1e18,'oc')
        py.plot(T_mir_arr, NEP[j,:]*1e18,'-c', label='140GHz')
        py.plot(T_mir_arr, NEP_ph_w[j,:]*1e18,'-.b')
        py.plot(T_mir_arr, NEP_th[j,:]*1e18,'--b')
        py.plot(T_mir_arr, NEP_readout[j,:]*1e18,'.b')

    if freq_str[j] == 195:
        py.plot(T_mir_arr, NEP[j,:]*1e18,'om')
        py.plot(T_mir_arr, NEP[j,:]*1e18,'-m', label='195GHz')
        py.plot(T_mir_arr, NEP_ph_w[j,:]*1e18,'-.b')
        py.plot(T_mir_arr, NEP_th[j,:]*1e18,'--b')
        py.plot(T_mir_arr, NEP_readout[j,:]*1e18,'.b')

    if freq_str[j] == 280:
        py.plot(T_mir_arr, NEP[j,:]*1e18,'oy')
        py.plot(T_mir_arr, NEP[j,:]*1e18,'-y', label='280GHz')
        py.plot(T_mir_arr, NEP_ph_w[j,:]*1e18,'-.b')
        py.plot(T_mir_arr, NEP_th[j,:]*1e18,'--b')
        py.plot(T_mir_arr, NEP_readout[j,:]*1e18,'.b')

py.xlabel('Temperature of optical system [K]')
py.ylabel('NEP [aW/rtHz]')
py.legend(loc='best')

py.show()
