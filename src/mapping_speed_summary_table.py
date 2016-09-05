import numpy as np
import pylab as py
import os 
import sys
import lib_mappingspeed as lib_ms
import global_par as g

pi = np.pi

#######################################################################################################
wafer_num = g.wafer_num
#;wafer_num = [1, 1, 1, 1, 1, 1]

#;if fix(total(wafer_num)) eq 6 then subdir_tag = '1wafer/'
#;if fix(total(wafer_num)) eq 39 then subdir_tag = 'HFP5_LFP8_wafers/'

aperture_diameter_mm = sys.argv[3] #; mm
freq_str = []
for i_str in g.freq_GHz_arr: freq_str.append(int(i_str))

num_freq = len(freq_str)

T_bath = g.T_bath
#T_mir_arr = [2., 2.5, 3.0, 3.5, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, \
#                 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.]
T_mir_arr = g.T_mir_arr
num_Tmir = len(T_mir_arr)

NEP = np.zeros((num_freq, num_Tmir))
NET = np.zeros((num_freq, num_Tmir))
Pload = np.zeros((num_freq, num_Tmir))
NEP_ph_w = np.zeros((num_freq, num_Tmir))
NEP_th = np.zeros((num_freq, num_Tmir))
NEP_readout = np.zeros((num_freq, num_Tmir))

dir_in = sys.argv[1]
dir_out = sys.argv[2]

for i in range(0,num_Tmir):
    dir = dir_in+'/'+str(int(T_bath*1e3))+'mK_Dapt'+str(int(aperture_diameter_mm))+'mm'
#    '/Users/tomotakematsumura/work/develop/LiteBIRD/NET/CalNET/results/20140717/'+str(int(T_bath*1e3)) \
#        +'mK_Dapt'+str(int(aperture_diameter_mm))+'mm'

    dir_data = dir+'/data'
    dir_out = dir+'/out'
    os.system('mkdir -p '+dir_out)
   
#; beam_ref = 2.481e-8 ; str
#; eff = 0.64
#;   restore, dir_data+'/eff.save'
#;   restore, dir_data+'/emiss.save'
#######################################################################################################
    sens_arr = []
    print ''
    print 'band  pixD apteff Tedge  Pload NEP_ph   NEP_th    NEP_readout  NEP_tot  NET     Wafer Npix/Wafer Ndet/Wafer Ntot NET_arr Sens.'
    print '[GHz] [mm] []    [dB]   [pW] [aW/rtHz] [aW/rtHz] [aW/rtHz]  [aW/rtHz] [uKrtsec] [#]   [#]        [#]        [#]  [uKrtsec] [uK.arcmin]'
    for j in range(0,num_freq):
        file_in = dir_data+'/mappingspeed_'+str(freq_str[j])+'GHz_F_T'+str(T_mir_arr[i])+'K.npy'
        output = np.load(file_in)
        # select the pixel size
        if freq_str[j] < g.freq_break: ind = np.where((output[10] > (g.pix_low*1e-3-0.0001)) & (output[10] < (g.pix_low*1e-3+0.0001)))
        if freq_str[j] > g.freq_break: ind = np.where((output[10] > (g.pix_high*1e-3-0.0001)) & (output[10] < (g.pix_high*1e-3+0.0001)))
        NEP_ph_w[j,i] = output[19,ind[0]]
        NEP_th[j,i] = output[20,ind[0]]
        NEP_readout[j,i] = output[37,ind[0]]
        NEP[j,i] = output[26,ind[0]]
        NET[j,i] = output[32,ind[0]]
        Pload[j,i] = output[17,ind[0]]
        apt_eff = output[5,ind[0]]
        Npix = output[9,ind[0]]
        if freq_str[j] > 110: Npix = 37
        Ndet_p_wafer = Npix*2.
        Ndet_tot = Npix*2.*wafer_num[j]

        fsky = g.fsky
        sens = (10800./pi * np.sqrt(8.*pi*fsky*NET[j,i]**2/(g.Tmis_sec*Ndet_tot)))
        sens_arr.append(sens)

        print freq_str[j], 
        print "%2.1f " % (output[10,ind[0]]*1e3), 
        print "%2.3f %2.1f " % (apt_eff, np.array(lib_ms.apteff2TedgeDB(apt_eff))),
        print "%2.3f" % (Pload[j,i]*1e12), 
        print "%2.3f" % (NEP_ph_w[j,i]*1e18), 
        print "%2.3f" % (NEP_th[j,i]*1e18), 
        print "%2.3f" % (NEP_readout[j,i]*1e18), 
        print "%2.3f" % (NEP[j,i]*1e18), 
        print "%2.3f" % (NET[j,i]),
        print "%2d" % (wafer_num[j]),
        print "%2d" % (Npix),
        print "%2d" % (Ndet_p_wafer),
        print "%2d" % (Ndet_tot),
        print "%2.3f" % (NET[j,i]/np.sqrt(float(Ndet_tot))),
        print "%2.2f" % sens

    print 'Total sensitivity %2.3f' % (1./np.sqrt(np.sum(1./np.array(sens_arr)**2)))
    print ''
    lib_ms.print_specifications(True)
sys.exit()

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
