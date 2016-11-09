import numpy as np
import pylab as py
import lib_mappingspeed as lib_ms
import global_par as g
import os
import sys

pi = np.pi
radeg = (180./pi)
c = 299792458.


#######################################################################################################
def mapping_speed( dir_out):

    print '######################################################################'
    print ' aperture diameter: ', g.aperture_diameter_mm
    print ' mirror1', g.T_mir1
    print ' mirror2', g.T_mir2
    print ' aperture', g.T_ape
    print ' hwp', g.T_hwp
    print ' 1K', g.T_1K
    print ' bath temperature:', g.T_bath
    print 'halfangle_edge_degs', g.halfangle_edge_degs
    print ''
    print 'design temperatures are'
    print ' mirror1', g.T_mir1_nominal
    print ' mirror2', g.T_mir2_nominal
    print ' aperture', g.T_ape_nominal
    print ' hwp', g.T_hwp_nominal
    print ' 1K', g.T_1K_nominal
    print ' bath', g.T_bath_nominal
    print '######################################################################'

  # define basic parameter
    
    dir_out = dir_out+'/'+str(int(g.T_bath*1e3))+'mK_Dapt'+str(int(g.aperture_diameter_mm))+'mm/data/'
    os.system('mkdir -p '+dir_out)
    print dir_out

    # CMB, HWP, Aperture, 2 Mirrors, 1K filter, Lens, Microstrip
    filename = dir_out+'emiss'
    np.save( filename+'.save', g.emiss_arr)

    # CMB, 
    # HWP(AR), 
    # Aperture (not include in here), 
    # 2 mirrors, 1K filter 0.95 (goal value, Ade filter is 0.9-0.95), 
    # Lens AR, microstrip and antenna eff from Toki's PB value
    filename = dir_out+'eff'
    np.save(filename+'.save', g.ref_arr)

    T_fp = {}; T_elements = {}
    T_fp['T_bath'] = g.T_bath
    T_fp['T_bath_nominal'] = g.T_bath_nominal

    T_elements['T_hwp'] = g.T_hwp
    T_elements['T_ape'] = g.T_ape
    T_elements['T_mir1'] = g.T_mir1
    T_elements['T_mir2'] = g.T_mir2
    T_elements['T_1K'] = g.T_1K
    T_elements['T_lenslet'] = g.T_lenslet

    T_elements['T_hwp_nominal'] = g.T_hwp_nominal
    T_elements['T_ape_nominal'] = g.T_ape_nominal
    T_elements['T_mir1_nominal'] = g.T_mir1_nominal
    T_elements['T_mir2_nominal'] = g.T_mir2_nominal
    T_elements['T_1K_nominal'] = g.T_1K_nominal
    T_elements['T_lenslet_nominal'] = g.T_lenslet_nominal

#    print 'e2e except aperture eff',  eff_e2e_eae 
#; define the basic parameters
    num = 1000

#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; LiteBIRD spec
    num_wafer = 1

    num_lensD = 400
    D_lens = (np.arange(num_lensD)+1)/float(num_lensD)*50.e-3 #; max 35mm
    Npix = lib_ms.pixel_count( D_lens, g.wafer_size, num_wafer)
    #Npix = [57.,57.,57.,57.,57.,57., 148.,111.,148.,111.,148.,111]
#    output = np.zeros((40,ntype))
    output_arr = {}

    Pmax = 0.
    num_band = len(g.freq_GHz_arr)
    output_arr = {}
    band_info = {}
    filename_arr = []
    NEP_tot_w = []
   
#    eff_arr = np.copy(g.ref_arr)
    eff_arr = np.copy(g.eff_arr)
    emiss_arr = np.copy(g.emiss_arr)
    ref_arr = np.copy(g.ref_arr)
#    print eff_arr
    for i in range(0, num_band):
        band_info['freq_GHz'] = g.freq_GHz_arr[i]
        band_info['bandwidth'] = g.bandwidth_arr[i]
       # band_info['NEP_tot_w'] = g.NEP_tot_arr[i]
        NEP_tot_w = g.NEP_tot_arr[i]

       # print 'NEP_tot_w=',NEP_tot_w[i]
        
        output = {}
        num_key_pre = 43
        output_d = np.zeros((num_key_pre,num_lensD))

        for j in range(0, num_lensD):
            num_pix = Npix[j]
            d_pixel_mm = D_lens[j]*1.e3

            freq_c = band_info['freq_GHz'] * 1.e9
            if g.option_proposal == 'US_MO_LFT':
                eff_arr[0] = 1. - g.ref_arr[0]
                eff_arr[1] = 1. - g.ref_arr[1] - (1. - np.exp(-2.*pi*freq_c/c*g.d_lens*g.n_hwp*g.losstan_hwp))
                eff_arr[2] = 1. - g.ref_arr[2]
                eff_arr[3] = 1. - g.ref_arr[3] - g.mirror_abs
                eff_arr[4] = 1. - g.ref_arr[4] - g.mirror_abs
                eff_arr[5] = 1. - g.ref_arr[5] - (1. - np.exp(-2.*pi*freq_c/c*g.d_filter*g.n_filter*g.losstan_filter))
                eff_arr[6] = 1. - g.ref_arr[6]
                eff_arr[7] = 1. - g.ref_arr[7]
                emiss_arr[0] = g.emiss_arr[0]
                emiss_arr[1] = 1. - np.exp(-2.*pi*freq_c/c*g.d_lens*g.n_hwp*g.losstan_hwp)
                emiss_arr[2] = g.emiss_arr[2]
                emiss_arr[3] = g.mirror_abs
                emiss_arr[4] = g.mirror_abs
                emiss_arr[5] = 1. - np.exp(-2.*pi*freq_c/c*g.d_filter*g.n_filter*g.losstan_filter)
                emiss_arr[6] = g.emiss_arr[6]
                emiss_arr[7] = g.emiss_arr[7]

            if g.option_proposal == 'US_CSR_LFT':
                eff_arr[0] = 1. - g.ref_arr[0]
                eff_arr[1] = 1. - g.ref_arr[1] - (1. - np.exp(-2.*pi*freq_c/c*g.d_lens*g.n_hwp*g.losstan_hwp))
                eff_arr[2] = 1. - g.ref_arr[2]
                eff_arr[3] = 1. - g.ref_arr[3] - g.mirror_abs
                eff_arr[4] = 1. - g.ref_arr[4] - g.mirror_abs
                eff_arr[5] = 1. - g.ref_arr[5] - (1. - np.exp(-2.*pi*freq_c/c*g.d_filter*g.n_filter*g.losstan_filter))
                eff_arr[6] = 1. - g.ref_arr[6]
                eff_arr[7] = 1. - g.ref_arr[7]
                emiss_arr[0] = g.emiss_arr[0]
                emiss_arr[1] = 1. - np.exp(-2.*pi*freq_c/c*g.d_lens*g.n_hwp*g.losstan_hwp)
                emiss_arr[2] = g.emiss_arr[2]
                emiss_arr[3] = g.mirror_abs
                emiss_arr[4] = g.mirror_abs
                emiss_arr[5] = 1. - np.exp(-2.*pi*freq_c/c*g.d_filter*g.n_filter*g.losstan_filter)
                emiss_arr[6] = g.emiss_arr[6]
                emiss_arr[7] = g.emiss_arr[7]


            if g.option_proposal == 'US_MO_HFT':
                eff_arr[0] = 1. - g.ref_arr[0]
                eff_arr[1] = 1. - g.ref_arr[1] - (1. - np.exp(-2.*pi*freq_c/c*g.d_lens*g.n_hwp*g.losstan_hwp))
                eff_arr[2] = 1. - g.ref_arr[2]
                eff_arr[3] = 1. - g.ref_arr[3] - (1. - np.exp(-2.*pi*freq_c/c*g.d_lens*g.n_lens*g.losstan_lens))
                eff_arr[4] = 1. - g.ref_arr[4] - (1. - np.exp(-2.*pi*freq_c/c*g.d_lens*g.n_lens*g.losstan_lens))
                eff_arr[5] = 1. - g.ref_arr[5] - (1. - np.exp(-2.*pi*freq_c/c*g.d_filter*g.n_filter*g.losstan_filter))
                eff_arr[6] = 1. - g.ref_arr[6]
                eff_arr[7] = 1. - g.ref_arr[7]
                emiss_arr[0] = g.emiss_arr[0]
                emiss_arr[1] = 1. - np.exp(-2.*pi*freq_c/c*g.d_lens*g.n_hwp*g.losstan_hwp)
                emiss_arr[2] = g.emiss_arr[2]
                emiss_arr[3] = 1. - np.exp(-2.*pi*freq_c/c*g.d_lens*g.n_lens*g.losstan_lens)
                emiss_arr[4] = 1. - np.exp(-2.*pi*freq_c/c*g.d_lens*g.n_lens*g.losstan_lens)
                emiss_arr[5] = 1. - np.exp(-2.*pi*freq_c/c*g.d_filter*g.n_filter*g.losstan_filter)
                emiss_arr[6] = g.emiss_arr[6]
                emiss_arr[7] = g.emiss_arr[7]

         #   print 'system efficiency except aperture eff'
#            num_elements = len(g.eff_arr)
#            eff_e2e_eae = 1.
#            for ii in range(0, num_elements):
#                eff_e2e_eae *= g.eff_arr[ii]

#            print eff_arr
#            sys.exit()

#            aperture_diameter_mm = g.aperture_diameter_mm
#            halfangle_edge_degs = g.halfangle_edge_degs
#            emiss_arr = g.emiss_arr
#            eff_arr = g.eff_arr
#            print eff_arr
#            print g.eff_arr
#            print ''
            output_per_dpix = lib_ms.mapping_speed_dfix( d_pixel_mm, Pmax, band_info, num_pix, \
                                                        g.aperture_diameter_mm, emiss_arr, eff_arr, \
                                                        T_fp, T_elements, \
                                                         g.halfangle_edge_degs )
 #           print g.eff_arr
#            sys.exit()

            num_key = len(output_per_dpix.keys())
            if num_key != num_key_pre: 
                print num_key, num_key_pre
                print 'the assumed number of keys is not the same as expected. Check mapping_speed_dfix!'
            keyname_arr = []
            output_key = np.zeros((num_key))
            for jj in range(num_key):
                keyname = output_per_dpix.keys()[jj]
                keyname_arr.append(keyname)
#                output_key[jj] = output_per_dpix[keyname]
                output_d[jj,j] = output_per_dpix[keyname]

        output = {'output': output_d, \
                  'keyname_arr': keyname_arr, \
                    'freq_GHz': g.freq_GHz_arr[i], \
                    'emiss_arr': emiss_arr, \
                    'eff_arr': eff_arr, \
                    'help':'output, freq_GHz'}
 
        filename = dir_out+'mappingspeed_'+str(int(g.freq_GHz_arr[i]))+'GHz'
        print filename
        np.save(filename, output)

        filename_arr.append(filename)
#        output_arr = {'freq_GHz':g.freq_GHz_arr, \
#                    'output_'+str(g.freq_GHz_arr[i]):output, \
#                    'help':'freq_GHz,output_#'}
    return filename_arr


