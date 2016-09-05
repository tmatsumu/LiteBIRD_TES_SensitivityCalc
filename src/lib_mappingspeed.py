import numpy as np
import pylab as py
import global_par as g
import os
import sys

pi = np.pi
radeg = (180./pi)
c = 299792458.

#;@/home/tmatsumu/codes/IDL/hwpsim/planckfunc.pro

#;#########################################################################
#; compute the mapping speed 
#; 2010-10-04: add the readout noise contribution based on 7pA/rt(Hz)*Vbias
#; 2013-6-23: change the thermal conductivity as n=3 for 300mK or higher
#; 2015-12-23: Found a bug in array assignment in mapping_speed 
#; 2015-12-23: revise significantly 

def print_specifications(option=True):
    if option == True:
        print '##########################################'
        print ''
        print 'LB specifications'
        print 'band center [GHz]', g.freq_GHz_arr
        print 'bandwidth', g.bandwidth
        print 'wafer size', g.wafer_size
        print 'emissivity per element', g.emiss
        print 'efficiency', g.eff
        print 'observational time [year]', g.Tmis_year
        print 'observational time [sec]', g.Tmis_sec
        print 'fractional sky', g.fsky
        print ''
        print 'detector specifications'
        print 'bath temperature', g.T_bath
        print 'wafer_num', g.wafer_num
        print 'current noise', g.current_noise
        print 'safety factor', g.X_fac
        print 'power low index of the thermal conductivity', g.n_tc
        print 'optimal ratio of bolo Tc to bath temperature', g.optimal_bolotemp_ratio
        print ''
        print 'optics specifications'
        print 'half-angle edge [degs]', g.halfangle_edge_degs
        print ''
        print 'Mirror temperature', g.T_mir_arr
        print 'temperature: aperture', g.T_ape
        print 'temperature: hwp', g.T_hwp
        print 'temperature: 1K', g.T_1K
        print ''
        print 'designed temperature: mirror', g.T_mir_nominal
        print 'designed temperature: aperture', g.T_ape_nominal
        print 'designed temperature: hwp', g.T_hwp_nominal
        print 'designed temperature: 1K', g.T_1K_nominal
        print ''
        print 'Jupiter'
        print 'beam_ref', g.beam_ref
        print 'planet temperature', g.T_p
        print 'reference beam size', g.ref_exp_beam
        print '##########################################'
        print ''

def function( Apt_eff, n):
#; this aperture efficiency is the
#; fraction of the power going through
#; the aperture given 
#; the n F lambda pixel size equivalent beam
#; given n=1, then apt_eff=0.45 which agrees to 
#; Adrian's 30GHz spill over.#
#
#;  y = 1.-exp(- alog(2.)*n^2./(1.028)^2. ) ; tomo original
#;  y = 1.-exp(- 0.6559*n^2. ) ; same as above, but in different format
    y = 1.-np.exp(- 0.73*n**2. ) #; based on Toki's equation
    return y


def Apt_eff_20130209( nu_GHz, D_aperture_mm, d_pix_mm):
    y = 1.-np.exp(-6.70e-12 * nu_GHz**2 * D_aperture_mm**2 * d_pix_mm**2)
    return y

def Apt_eff_20140917( Tedge):
    return 1.-np.exp(-2.*0.115*Tedge)

def aperture(dpix_m,nu_Hz,halfangle_edge_deg):
#    sigma0_rad = 1./2.*g.c/nu/pi*2.6/dpix
    sigma0_rad = 1./2.*g.c/(nu_Hz)/pi*2.6/(dpix_m)
    Tedge_dB = - 10.* np.log10( np.exp(-(halfangle_edge_deg/radeg)**2/2. * (1./sigma0_rad)**2) )
    Pap = 1.-np.exp(-2.*0.115*Tedge_dB)
#    print ''
#    print 'aperture'
#    print g.c, nu_Hz*1e-9, pi, dpix_m*1e3, sigma0_rad, sigma0_rad*radeg
#    print halfangle_edge_deg, halfangle_edge_deg/radeg
#    print Tedge_dB
#    sys.exit()
    return sigma0_rad, Tedge_dB, Pap

def aperture_beamwaistfactor(dpix_m,nu_Hz,halfangle_edge_deg,beamwaistfactor):
#    sigma0_rad = 1./2.*g.c/nu/pi*2.6/dpix
    sigma0_rad = 1./2.*g.c/(nu_Hz)/pi*beamwaistfactor/(dpix_m)
    Tedge_dB = - 10.* np.log10( np.exp(-(halfangle_edge_deg/radeg)**2/2. * (1./sigma0_rad)**2) )
    constant = np.log(10.)/10.
    Pap = 1.-np.exp(-constant*Tedge_dB)
#    print ''
#    print 'aperture'
#    print g.c, nu_Hz*1e-9, pi, dpix_m*1e3, sigma0_rad, sigma0_rad*radeg
#    print halfangle_edge_deg, halfangle_edge_deg/radeg
#    print Tedge_dB
#    sys.exit()
    return sigma0_rad, Tedge_dB, Pap

def apteff2TedgeDB(apt_eff):
    return np.log(1.-apt_eff)/(-2.*0.115)

def Brightness_BB( f, T): #; unit of [W/m^2/sr/Hz]
    out = 2. * (g.h*f**3.)/g.c**2. /(np.exp( g.h/g.k_b*f/T) - 1.)
    return out

def flux( freq, T):
    out =  g.h*freq /(np.exp( g.h/g.k_b*freq/T) - 1.)
    return out

def Brightness_dBdT( freq, T): #; unit of [W/m^2/sr/Hz]
    x = g.h/g.k_b*freq/T
    prefac = 2./(T**2.*g.c*2.)*( (g.h*freq**2.)**2/g.k_b )
    dBdT = prefac * (-np.exp(x)/(np.exp(x)-1.)**2.)
    return dBdT

#;-----
def Brightness_dPdT( freq_GHz, T, AOmega, eff, bandwidth): #; unit of [W/m**2/sr/Hz]
    num = 1000
    freq_i = freq_GHz*(1.-bandwidth/2.)*1.e9
    freq_f = freq_GHz*(1.+bandwidth/2.)*1.e9 
    freq_arr = freq_i+np.arange(num)/float(num)*(freq_f-freq_i)
    x = g.h/g.k_b*freq_arr/T
    prefac = 2./(T**2.*g.c**2.)*( (g.h*freq_arr**2.)**2/g.k_b )
    dBdT = 0.5 * AOmega * prefac * (-np.exp(x)/(np.exp(x)-1.)**2.) *eff
    int_dPdT = np.sum( dBdT*(freq_arr[num-1]-freq_arr[0])/float(num))
    return int_dPdT


def NEP2NET_CMB( NEP, freq, AOmega, eff):
    num = len(freq)
    x = g.h/g.k_b*freq/g.Tcmb
    prefac = 2./(g.Tcmb**2.*g.c**2.)*( (g.h*freq**2.)**2/g.k_b )
    dPdT = 0.5 * AOmega * prefac * (-np.exp(x)/(np.exp(x)-1.)**2.) *eff
    int_dPdT = np.sum( dPdT * (freq[num-1]-freq[0]) / float(num) )
    #  int_dPdT = int_tabulated( freq, dPdT, /double) 
#;  print, '->', total( 0.5 * AOmega * prefac * (-np.exp(x)/(np.exp(x)-1.)**2.) * (freq[num-1]-freq[0]) / float(num) )
    NET2 = 1. / int_dPdT**2. * NEP**2.
    return np.sqrt(NET2)/np.sqrt(2.)

def Planet_dilution_T2P( freq_GHz, T_p, beam_ref, beam_exp, bandwidth):
    #; beam_ref, beam_exp in [str]
    num = 1000
    freq_i = freq_GHz*(1.-bandwidth/2.)*1.e9
    freq_f = freq_GHz*(1.+bandwidth/2.)*1.e9 
    freq = freq_i+np.arange(num)/float(num)*(freq_f-freq_i)
    A_omega = (g.c/freq)**2.
    BB = Brightness_BB(freq, T_p)                                #; 1.86K at 143GHz
#    P = 0.5*int_tabulated( freq, A_omega*BB, /double)    #;*eff_opt
    P = 0.5*np.sum(A_omega*BB*(freq[num-1]-freq[0])/float(num))
    P = P*beam_ref/beam_exp
    return P

def pixel_count( pixel_diam, s2s_wafer_size, num_wafer):
    Npix = np.int_(num_wafer*np.sqrt(3.)/2.*s2s_wafer_size**2./6.*pi*np.sqrt(3.)*4./pi/pixel_diam**2)
    return Npix

#;####################################################################################################################################


def mapping_speed_dfix(d_pixel_mm, Pmax, band_info, Npix, aperture_diameter_mm, emiss_arr, eff_arr, \
                       T_fp, T_elements, \
                       halfangle_edge_degs):
    '''
    d_pixel_mm, Pmax, band_info, Npix, aperture_diameter_mm, emiss_arr, eff_arr, \
                           T_fp, T_elements, \
                           halfangle_edge_degs
    '''
    freq_GHz = band_info['freq_GHz']
    bandwidth = band_info['bandwidth']

    T_bath = T_fp['T_bath']
    T_bath_nominal = T_fp['T_bath_nominal']

    T_hwp = T_elements['T_hwp']
    T_ape = T_elements['T_ape']
    T_mir1 = T_elements['T_mir1']
    T_mir2 = T_elements['T_mir2']
    T_1K = T_elements['T_1K']
    T_lenslet = T_elements['T_lenslet']


    T_hwp_nominal = T_elements['T_hwp_nominal']
    T_ape_nominal = T_elements['T_ape_nominal']
    T_mir1_nominal = T_elements['T_mir1_nominal']
    T_mir2_nominal = T_elements['T_mir2_nominal']
    T_1K_nominal = T_elements['T_1K_nominal']
    T_lenslet_nominal = T_elements['T_lenslet_nominal']

    num_elements = len(eff_arr)
    eff_loc = np.zeros(num_elements)
    eff_e2e_eae = 1.
#    print 'system efficiency except aperture eff'
    for i in range(0, num_elements):
        eff_e2e_eae *= eff_arr[i]
        eff_loc[i] = eff_arr[i]
#        if d_pixel_mm == 18: print '->', eff_loc[i], eff_arr[i], eff_e2e_eae
#    print 'e2e except aperture eff',  eff_e2e_eae 
#; define the basic parameters
    num = 1000

#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; LiteBIRD spec

    T_fil = g.T_1K_nominal
    num_wafer = 1
    D_lens = d_pixel_mm*1.e-3 # (np.arange(400)+1)/400.*50.e-3 #; max 35mm
    freq_c = float(freq_GHz)*1.e9
#    lambda = c/freq_c
    output = np.zeros(40)
#    sigma_0, Tedge_dB, apt_eff = aperture(D_lens,freq_c,halfangle_edge_degs)
    sigma_0, Tedge_dB, apt_eff = aperture_beamwaistfactor(D_lens,freq_c,halfangle_edge_degs,g.beamwaistfactor)

    eff_e2e = eff_e2e_eae*apt_eff
    num_bolo = 2.*Npix
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; OPTICAL LOADING AND NEP_optical
 
#; define the frequency range
    freq_i = freq_c*(1.-bandwidth/2.)
    freq_f = freq_c*(1.+bandwidth/2.)
    freq = freq_i + np.arange(num)/float(num)*(freq_f-freq_i)
    A_omega = (g.c/freq)**2.

#; integrate over the band
#;  CMB only
    BB = Brightness_BB(freq, g.Tcmb)
    F_cmb = flux(freq,g.Tcmb)
    P_cmb = np.sum( A_omega*BB * (freq[num-1]-freq[0])/float(num)) # ;*eff_opt*eff_spill*T_trans
    
#; HWP ony
    BB = Brightness_BB(freq, T_hwp)
    F_hwp = flux(freq,T_hwp)
    P_hwp = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*(1.-eff_spill)*T_trans

#; Nominal HWP temperature
    BB = Brightness_BB(freq, T_hwp_nominal)
    F_hwp_nominal = flux(freq,T_hwp_nominal)
    P_hwp_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #*eff_opt*(1.-eff_spill)*T_trans
      
#; Aperture only
    BB = Brightness_BB(freq, T_ape)
    F_aperture = flux(freq,T_ape)
    P_aperture = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*(1.-eff_spill)*T_trans

#; Nominal Aperture temperature
    BB = Brightness_BB(freq, T_ape_nominal)
    F_aperture_nominal = flux(freq,T_ape_nominal)
    P_aperture_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #*eff_opt*(1.-eff_spill)*T_trans
      
#; Mirror1 only
    BB = Brightness_BB(freq, T_mir1)
    F_mir1 = flux(freq,T_mir1)
    P_mir1 = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) # ;*eff_opt*emiss_mirror*2.*T_trans

#; Nominal Mirror1 temperature 
    BB = Brightness_BB(freq, T_mir1_nominal)
    F_mir1_nominal = flux(freq,T_mir1_nominal)
    P_mir1_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans

#; Mirror2 only
    BB = Brightness_BB(freq, T_mir2)
    F_mir2 = flux(freq,T_mir2)
    P_mir2 = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) # ;*eff_opt*emiss_mirror*2.*T_trans

#; Nominal Mirror2 temperature 
    BB = Brightness_BB(freq, T_mir2_nominal)
    F_mir2_nominal = flux(freq,T_mir2_nominal)
    P_mir2_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans

#; 1K filter only
    BB = Brightness_BB(freq, T_fil)
    F_1Kfilter = flux(freq,T_fil)
    P_1Kfilter = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans

#; nominal 1K filter temperature
    BB = Brightness_BB(freq, T_1K_nominal)
    F_1Kfilter_nominal = flux(freq,T_1K_nominal)
    P_1Kfilter_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans

#; lenslet temperature
    BB = Brightness_BB(freq, T_lenslet)
    F_lenslet = flux(freq,T_lenslet)
    P_lenslet = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans

#; nominal lenslet temperature
    BB = Brightness_BB(freq, T_lenslet_nominal)
    F_lenslet_nominal = flux(freq,T_lenslet_nominal)
    P_lenslet_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans
 
#; Pixel lens only
    BB = Brightness_BB(freq, T_bath)
    F_100mK = flux(freq, T_bath)
    P_100mK =  np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt

#; nominal Pixel lens temperature
    BB = Brightness_BB(freq, T_bath_nominal)
    F_100mK_nominal = flux(freq, T_bath_nominal)
    P_100mK_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt
 
#; Jupiter
    BB = Brightness_BB(freq, 1.86)                          #; 1.86K at 143GHz
    F_Jupiter = flux(freq,1.86)
    P_jupiter = 0.5*np.sum(A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt

    eff_loc[2] = apt_eff

    F_sum =   ( F_cmb     * emiss_arr[0] *eff_loc[1]*eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  F_hwp      * emiss_arr[1]            *eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  F_aperture * emiss_arr[2] *(1.-eff_loc[2])      *eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  F_mir1     * emiss_arr[3]                                  *eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  F_mir2     * emiss_arr[4]                                             *eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  F_1Kfilter * emiss_arr[5]                                                        *eff_loc[6]*eff_loc[7] \
            +  F_lenslet  * emiss_arr[6]                                                                   *eff_loc[7] )

    F2_sum =  ( F_cmb      * emiss_arr[0] *eff_loc[1]*eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] )**2 \
            + ( F_hwp      * emiss_arr[1]            *eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] )**2 \
            + ( F_aperture * emiss_arr[2] *(1.-eff_loc[2])      *eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] )**2 \
            + ( F_mir1     * emiss_arr[3]                                  *eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] )**2 \
            + ( F_mir2     * emiss_arr[4]                                             *eff_loc[5]*eff_loc[6]*eff_loc[7] )**2 \
            + ( F_1Kfilter * emiss_arr[5]                                                        *eff_loc[6]*eff_loc[7] )**2 \
            + ( F_lenslet  * emiss_arr[6]                                                                   *eff_loc[7] )**2


    P_load = 0.5 \
            * ( P_cmb     * emiss_arr[0] *eff_loc[1]*eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  P_hwp      * emiss_arr[1]            *eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  P_aperture * emiss_arr[2] *(1.-eff_loc[2])      *eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  P_mir1     * emiss_arr[3]                                  *eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  P_mir2     * emiss_arr[4]                                             *eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  P_1Kfilter * emiss_arr[5]                                                        *eff_loc[6]*eff_loc[7] \
            +  P_lenslet  * emiss_arr[6]                                                        *eff_loc[6]*eff_loc[7] )

    P_load_nominal = 0.5 \
            * ( P_cmb     * emiss_arr[0] *eff_loc[1]*eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  P_hwp      * emiss_arr[1]            *eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  P_aperture * emiss_arr[2] *(1.-eff_loc[2])      *eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  P_mir1     * emiss_arr[3]                                  *eff_loc[4]*eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  P_mir2     * emiss_arr[4]                                             *eff_loc[5]*eff_loc[6]*eff_loc[7] \
            +  P_1Kfilter * emiss_arr[5]                                                        *eff_loc[6]*eff_loc[7] \
            +  P_lenslet  * emiss_arr[6]                                                        *eff_loc[6]*eff_loc[7] )

    # F_sum = eff_loc[6]*eff_loc[5] \
    #      *(F_cmb*apt_eff_arr*eff_loc[0]*emiss[0] *eff_loc[1]*eff_loc[3]*eff_loc[4] \
    #        +F_hwp*apt_eff_arr*eff_loc[1]*emiss[1] *eff_loc[3]*eff_loc[4] \
    #        +F_aperture*(1.-apt_eff_arr)*emiss[2] *eff_loc[3]*eff_loc[4] \
    #        +F_mir1*2.*eff_loc[3]*emiss[3]  *emiss[4]\
    #        +F_1Kfilter*eff_loc[4]*emiss[4] )

    # F2_sum = (eff_loc[6]*eff_loc[5])**2  \
    #      * ( (F_cmb*apt_eff_arr*eff_loc[0]*emiss[0]  *eff_loc[1]*eff_loc[3]*eff_loc[4])**2  \
    #              +(F_hwp*apt_eff_arr*eff_loc[1]*emiss[1] *eff_loc[3]*eff_loc[4])**2  \
    #              +(F_aperture*(1.-apt_eff_arr)*emiss[2]   *eff_loc[3]*eff_loc[4])**2  \
    #              +(F_mir1*2.*eff_loc[3]*emiss[3]  *eff_loc[4])**2  \
    #              +(F_1Kfilter*eff_loc[4]*emiss[4])**2 )              
    
#    P_load = 0.5*eff_loc[6]*eff_loc[5]  \
#         * (P_cmb*apt_eff_arr*eff_loc[0]*emiss[0] *eff_loc[1]*eff_loc[3]*eff_loc[4] \
#                + P_hwp*apt_eff_arr*eff_loc[1]*emiss[1]  *eff_loc[3]*eff_loc[4]\
#                + P_aperture*(1.-apt_eff_arr)*emiss[2]  *eff_loc[3]*eff_loc[4]\
#                + P_mir1*2.*eff_loc[3]*emiss[3]  *eff_loc[4]\
#                + P_1Kfilter*eff_loc[4]*emiss[4] )

    # P_load_nominal = 0.5*eff_loc[6]*eff_loc[5]  \
    #      * (P_cmb*apt_eff_arr*eff_loc[0]*emiss[0]   *eff_loc[1]*eff_loc[3]*eff_loc[4]\
    #             + P_hwp_nominal*apt_eff_arr*eff_loc[1]*emiss[1] *eff_loc[3]*eff_loc[4] \
    #             + P_aperture_nominal*(1.-apt_eff_arr)*emiss[2]   *eff_loc[3]*eff_loc[4]\
    #             + P_mir1_nominal*2.*eff_loc[3]*emiss[3]  *eff_loc[4] \
    #             + P_1Kfilter_nominal*eff_loc[4]*emiss[4] )

    poisson2 = np.sum( 2.*F_sum * g.h * freq * (freq[num-1]-freq[0])/float(num))
#        bunching2 = np.sum(F2_sum * (freq[num-1]-freq[0])/float(num))
#        bunching2 = np.sum(2.*F2_sum * (freq[num-1]-freq[0])/float(num))
    bunching2 = np.sum(2.*F_sum**2. * (freq[num-1]-freq[0])/float(num))
 
    NEP_ph_w = np.sqrt(poisson2 + bunching2)
    NEP_ph_wo = np.sqrt(poisson2)
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; THERMAL NOISE

    if (P_load > Pmax): P_load_max = P_load_nominal
    if (P_load <= Pmax): 
        P_load_max = Pmax
        print 'maximum loading is taken for Jupiter', \
        freq_c*1e-9, '%2.3f %2.3f %2.3f' % (P_load_nominal*1e12, P_load_max*1e12, Pmax*1e12)
    P_load_max = P_load_nominal
    n_tc = g.n_tc

    optimal_bolotemp_ratio = g.optimal_bolotemp_ratio
    T_bolo = optimal_bolotemp_ratio*T_bath_nominal      #    ;0.2
    X_fac = g.X_fac  # Psat fact 3

    gamma = (n_tc+1.)/(2.*n_tc+3.)*(1.-(T_bath/T_bolo)**(2.*n_tc+3.))/(1.-(T_bath/T_bolo)**(n_tc+1.))
    Gave = X_fac*P_load_max/(T_bolo-T_bath)
    NEP_th_2 = 4.*g.k_b*X_fac*P_load_max*T_bath   \
        * (n_tc+1.)**2/(2.*n_tc+3.)  \
        * ((T_bolo/T_bath)**(2.*n_tc+3)-1.) / ((T_bolo/T_bath)**(n_tc+1)-1.)**2 

    NEP_th = np.sqrt(NEP_th_2)
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; READOUT NOISE
    R_bolo = 1.
    Vbias = np.sqrt((X_fac-1.)*R_bolo*P_load_max)
    NEP_readout = np.sqrt(NEP_ph_w**2. + NEP_th**2.)*np.sqrt(1.1**2-1.)
#    NEP_readout = g.current_noise*Vbias

#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; NEP_total     
#    NEP_to_w = np.sqrt(NEP_ph_w**2. + NEP_th**2. + NEP_readout**2.)    #  ;*np.sqrt(2.)
#    NEP_to_wo = np.sqrt(NEP_ph_wo**2. + NEP_th**2. + NEP_readout**2.)  #  ;*np.sqrt(2.)
    NEP_to_w = np.sqrt(NEP_ph_w**2. + NEP_th**2.)*1.1    #  ;*np.sqrt(2.)
    NEP_to_wo = np.sqrt(NEP_ph_wo**2. + NEP_th**2.)*1.1  #  ;*np.sqrt(2.)

#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; SUMMARY and PRINT
#     
#; convert from NEP to NET
    NET_photon_w= NEP2NET_CMB(NEP_ph_w,freq,A_omega,eff_e2e)*1.e6
    NET_photon_wo= NEP2NET_CMB(NEP_ph_wo,freq,A_omega,eff_e2e)*1.e6
    NET_detect= NEP2NET_CMB(NEP_th,freq,A_omega,eff_e2e)*1.e6
    NET_readout= NEP2NET_CMB(NEP_readout,freq,A_omega,eff_e2e)*1.e6
    NET_totalper_w = NEP2NET_CMB(NEP_to_w,freq,A_omega,eff_e2e)*1.e6
    NET_totalper_wo = NEP2NET_CMB(NEP_to_wo,freq,A_omega,eff_e2e)*1.e6
    NEQ_total_w = (np.sqrt(2.)*NET_totalper_w)/np.sqrt(float(Npix))
    NEQ_total_wo = (np.sqrt(2.)*NET_totalper_wo)/np.sqrt(float(Npix))

#     if d_pixel_mm == 18:
#         print freq_GHz
#         print ''
#         print P_cmb*1e12, P_hwp*1e12, P_aperture*1e12, P_mir1*1e12, P_mir2*1e12, P_1Kfilter*1e12
#         print 'P_cmb', (P_cmb     * emiss_arr[0] *eff_loc[1]*eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6])*1e12
#         print emiss_arr[0] *eff_loc[1]*eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]
#         print emiss_arr[0]
#         print eff_loc[1]*eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6]
#         print eff_loc[1:6]
#         print apt_eff
#         print ''
#         print 'P_hwp', (P_hwp      * emiss_arr[1]        *eff_loc[2]*eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6])*1e12
#         print 'P_aperture', (P_aperture * emiss_arr[2] *(1.-eff_loc[2])  *eff_loc[3]*eff_loc[4]*eff_loc[5]*eff_loc[6])*1e12
#         print 'P_mir1', (P_mir1     * emiss_arr[3]                      *eff_loc[4]*eff_loc[5]*eff_loc[6])*1e12
#         print 'P_mir2', (P_mir2     * emiss_arr[4]                             *eff_loc[5]*eff_loc[6])*1e12
#         print 'P_1Kfilter', (P_1Kfilter * emiss_arr[5]                                    *eff_loc[6])*1e12
#         print 'P_load', P_load*1e12
#         print ''
#         print NEP_to_w
#         print NEP2NET_CMB(NEP_to_w,freq,A_omega,eff_e2e)*1.e6
#         print NEP2NET_CMB(1.,freq,A_omega,1.)*1.e6
#         print NEP2NET_CMB(1.,freq,A_omega,eff_e2e)*1.e6
#         print np.mean(A_omega), np.mean(np.sqrt(A_omega))
#         print eff_e2e_eae
#         print eff_e2e
# #        sys.exit()

    output = {'freq_c':freq_c, \
            'freq_i': freq_i, \
            'freq_f': freq_f, \
            'bandwidth': bandwidth, \
            'mean_AOmega': np.mean(A_omega), \
            'apt_eff': apt_eff, \
            'T_mir1': T_mir1, \
            'T_mir2': T_mir2, \
            'Npix': Npix, \
            'D_lens': D_lens, \
            'Pmax': Pmax, \
            'P_cmb': P_cmb, \
            'P_hwp': P_hwp, \
            'P_aperture': P_aperture, \
            'P_mir1': P_mir1, \
            'P_mir2': P_mir2, 
            'P_1Kfilter': P_1Kfilter, 
            'P_100mK': P_100mK, \
            'P_jupier': P_jupiter, 
            'P_load': P_load, \
            'T_bath': T_bath, \
            'n_tc': n_tc, \
            'T_bolo': T_bolo, \
            'gamma': gamma, \
            'X_fac': X_fac, \
            'Gave': Gave, \
            'Vbias': Vbias, \
            'NEP_ph_wo': NEP_ph_wo, \
            'NEP_ph_w': NEP_ph_w, \
            'NEP_th': NEP_th, \
            'NEP_readout': NEP_readout, \
            'NEP_to_w': NEP_to_w, \
            'NEP_to_wo': NEP_to_wo, \
            'NET_photon_w': NET_photon_w, \
            'NET_photon_wo': NET_photon_wo, \
            'NET_detect': NET_detect, \
            'NET_readout': NET_readout, \
            'NET_totalper_w': NET_totalper_w, \
            'NET_totalper_wo': NET_totalper_wo, \
            'NEQ_total_w': NEQ_total_w, \
            'NEQ_total_wo': NEQ_total_wo, \
            'NEP2NET_CMB': NEP2NET_CMB(1.,freq,A_omega,eff_e2e)*1.e6, \
            'eff_e2e_eae': eff_e2e_eae} 
#            'eff_loc': eff_loc, \
#            'emiss_arr': emiss_arr }
#            'help': 'freq_c,freq_i,freq_f,bandwidth,mean_AOmega,apt_eff,T_mir1,T_mir2,Npix,D_lens,Pmax,P_cmb,P_hwp,P_aperture,P_mir1,P_mir2,P_1Kfilter,P_100mK,P_jupiter,P_load,T_bath,n_tc,T_bolo,gamma,X_fac,Gave,NEP_ph_wo,NEP_ph_w,NEP_th,NEP_readout,NEP_to_w,NEP_to_wo,NET_photon_w,NET_photon_wo,NET_detect,NET_readout,NET_totalper_w,NET_totalper_wo,NEQ_total_w,NEQ_total_wo,NEP2NET_CMB,eff_e2e_eae'}
    return output
#end

#;####################################################################################################################################
