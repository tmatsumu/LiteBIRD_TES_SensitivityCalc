import numpy as np
import pylab as py
import global_par as g
import os
import sys

pi = np.pi
radeg = (180./pi)

#;@/home/tmatsumu/codes/IDL/hwpsim/planckfunc.pro

#;#########################################################################
#; compute the mapping speed 
#; 2010-10-04: add the readout noise contribution based on 7pA/rt(Hz)*Vbias
#; 2013-6-23: change the thermal conductivity as n=3 for 300mK or higher

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

def mapping_speed( freq_GHz, Pmax, dir_out, bandwidth, wafer_size, aperture_diameter_mm, emiss, eff, \
                       T_bath, T_mir, T_ape, T_hwp, T_1K, \
                       T_bath_nominal, T_mir_nominal, T_ape_nominal, T_hwp_nominal, T_1K_nominal, halfangle_edge_degs):

    num_elements = len(eff)
    eff_e2e_eae = 1.
    print 'system efficiency except aperture eff'
    for i in range(0, num_elements):
        eff_e2e_eae *= eff[i]
    print 'e2e except aperture eff',  eff_e2e_eae 
#; define the basic parameters
    num = 1000

#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; LiteBIRD spec

    T_fil = g.T_1K_nominal
#;  T_mir_nominal = 4.
#;  T_ape_nominal = 4.
#;  T_HWP_nominal = 4.
#;  T_bath_nominal = 0.1
#;  T_1K_nominal = 1.

    num_wafer = 1

    D_lens = (np.arange(400)+1)/400.*50.e-3 #; max 35mm
    Npix = pixel_count( D_lens, wafer_size, num_wafer)
  
    freq_c = float(freq_GHz)*1.e9
#    lambda = c/freq_c
    ntype=len(Npix)

    output = np.zeros((40,ntype))

    n = (np.arange(400)+1)/400.*10. ; 

    for j in range(0, ntype):
        num_pix = Npix[j]
        d_pixel_mm = D_lens[j]*1.e3
#        apt_eff = Apt_eff_20130209(freq_GHz, aperture_diameter_mm, d_pixel_mm)
#        apt_eff = Apt_eff_20140916(Tedge)
#        apt_eff = Apt_eff_20131206(freq_GHz, n[j])

        sigma_0, Tedge_dB, apt_eff = aperture(d_pixel_mm*1.e-3,freq_c,halfangle_edge_degs)
#        sigma_0, Tedge_dB, apt_eff = g.aperture_ext(freq_c,'Toki')

        eff_e2e = eff_e2e_eae*apt_eff
#        print apt_eff, eff_e2e
#        sys.exit()
        
        num_bolo = 2.*Npix
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; OPTICAL LOADING AND NEP_optical
     
#; define the frequency range
        freq_i = freq_c*(1.-bandwidth/2.)
        freq_f = freq_c*(1.+bandwidth/2.)
        freq = freq_i + np.arange(num)/float(num)*(freq_f-freq_i)
        A_omega = (g.c/freq)**2.
        output[0,j]=freq_c
        output[1,j]=freq_i
        output[2,j]=freq_f
        output[3,j]=bandwidth
        output[4,j]=np.mean(A_omega)
        output[5,j]=apt_eff

        output[6,j] = wafer_size
        output[7,j] = T_mir
        output[8,j] = Pmax
        output[9,j] = Npix[j]
        output[10,j] = D_lens[j]

#; integrate over the band
#;  CMB only
        BB = Brightness_BB(freq, g.Tcmb)
        F_cmb = flux(freq,g.Tcmb)
        P_cmb = np.sum( A_omega*BB * (freq[num-1]-freq[0])/float(num)) # ;*eff_opt*eff_spill*T_trans
        output[11,j]=P_cmb
        
#; HWP ony
        BB = Brightness_BB(freq, T_hwp)
        F_hwp = flux(freq,T_hwp)
        P_hwp = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*(1.-eff_spill)*T_trans
        output[12,j]=P_hwp

#; Nominal HWP temperature
        BB = Brightness_BB(freq, T_hwp_nominal)
        F_hwp_nominal = flux(freq,T_hwp_nominal)
        P_hwp_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #*eff_opt*(1.-eff_spill)*T_trans
          
#; Aperture only
        BB = Brightness_BB(freq, T_ape)
        F_aperture = flux(freq,T_ape)
        P_aperture = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*(1.-eff_spill)*T_trans
        output[12,j]=P_aperture

#; Nominal Aperture temperature
        BB = Brightness_BB(freq, T_ape_nominal)
        F_aperture_nominal = flux(freq,T_ape_nominal)
        P_aperture_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #*eff_opt*(1.-eff_spill)*T_trans
          
#; Mirror only
        BB = Brightness_BB(freq, T_mir)
        F_mirror = flux(freq,T_mir)
        P_mirror = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) # ;*eff_opt*emiss_mirror*2.*T_trans
        output[13,j]=P_mirror

#; Nominal Mirror temperature 
        BB = Brightness_BB(freq, T_mir_nominal)
        F_mirror_nominal = flux(freq,T_mir_nominal)
        P_mirror_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans

#; 1K filter only
        BB = Brightness_BB(freq, T_fil)
        F_1Kfilter = flux(freq,T_fil)
        P_1Kfilter = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans
        output[14,j]=P_1Kfilter

#; nominal 1K filter temperature
        BB = Brightness_BB(freq, T_1K_nominal)
        F_1Kfilter_nominal = flux(freq,T_1K_nominal)
        P_1Kfilter_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans
     
#; Pixel lens only
        BB = Brightness_BB(freq, T_bath)
        F_100mK = flux(freq, T_bath)
        P_100mK =  np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt
        output[15,j]=P_100mK

#; nominal Pixel lens temperature
        BB = Brightness_BB(freq, T_bath_nominal)
        F_100mK_nominal = flux(freq, T_bath_nominal)
        P_100mK_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt
     
#; Jupiter
        BB = Brightness_BB(freq, 1.86)                          #; 1.86K at 143GHz
        F_Jupiter = flux(freq,1.86)
        P_jupiter = 0.5*np.sum(A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt
        output[16,j]=P_jupiter
     
#; SET THE HWP AND APERTURE THE SAME TEMPERATURE
#; Sum all the powers

# Microstrip related efficiency and lenslet AR
# CMB, aperture eff
# HWP, emissivity, aperture_eff, AR
# Aperture, 1-aperture_eff
# Mirrors, mirror_emiss
# 1K filter, emissivity, efficiency

#        print ''
#        print 'overall', eff[6]
#        print 'Lens let AR', eff[5]
#        print 'source  T[K],  apt_eff,  eff,  emis, P_each [pW], (P_each x apt_eff x eff x emiss) [pW]' 
#        print 'sky signal: %2.2f %2.2f %2.2f %2.2f %2.3f %2.3f' % (g.Tcmb, apt_eff, eff[0], emiss[0], P_cmb*1e12, 0.5*eff[6]*eff[5] * (P_cmb*apt_eff*eff[0]*emiss[0])*1e12 )
#        print 'HWP: %2.1f %2.2f %2.2f %2.2f %2.3f %2.3f' % (T_hwp, apt_eff, eff[1], emiss[1], P_hwp*1e12, 0.5*eff[6]*eff[5] * P_hwp*apt_eff*eff[1]*emiss[1]*1e12)
#        print 'Lyot: %2.1f %2.2f %2.2f %2.2f %2.3f %2.3f' % (T_ape, (1.-apt_eff), 1., emiss[2], P_aperture*1e12, 0.5*eff[6]*eff[5] * P_aperture*(1.-apt_eff)*emiss[2]*1e12)
#        print '2 mirrors: %2.1f %2.2f %2.2f %2.2f %2.3f %2.3f' % (T_mir, 1.,eff[3], emiss[3], P_mirror*1e12, 0.5*eff[6]*eff[5] * P_mirror*2.*eff[3]*emiss[3]*1e12)
#        print 'Filter: %2.1f %2.2f %2.2f %2.2f %2.3f %2.3f' % (T_fil, 1.,eff[4], emiss[4], P_1Kfilter*1e12, 0.5*eff[6]*eff[5] * P_1Kfilter*eff[4]*emiss[4]*1e12)
#        sys.exit()

        F_sum = eff[6]*eff[5] \
            *(F_cmb*apt_eff*eff[0]*emiss[0] *eff[1]*eff[3]*eff[4] \
              +F_hwp*apt_eff*eff[1]*emiss[1] *eff[3]*eff[4] \
              +F_aperture*(1.-apt_eff)*emiss[2] *eff[3]*eff[4] \
              +F_mirror*2.*eff[3]*emiss[3]  *emiss[4]\
              +F_1Kfilter*eff[4]*emiss[4] )

# Microstrip related efficiency and lenslet AR
# CMB, aperture eff
# HWP, emissivity, aperture_eff, AR
# Aperture, 1-aperture_eff
# Mirrors, mirror_emiss
# 1K filter, emissivity, efficiency
        F2_sum = (eff[6]*eff[5])**2  \
            * ( (F_cmb*apt_eff*eff[0]*emiss[0]  *eff[1]*eff[3]*eff[4])**2  \
                    +(F_hwp*apt_eff*eff[1]*emiss[1] *eff[3]*eff[4])**2  \
                    +(F_aperture*(1.-apt_eff)*emiss[2]   *eff[3]*eff[4])**2  \
                    +(F_mirror*2.*eff[3]*emiss[3]  *eff[4])**2  \
                    +(F_1Kfilter*eff[4]*emiss[4])**2 )              
        
# Microstrip related efficiency and lenslet AR
# CMB, aperture
# HWP, emissivity, apt_eff, AR
# Aperture, 1-apt_eff
        P_load = 0.5*eff[6]*eff[5]  \
            * (P_cmb*apt_eff*eff[0]*emiss[0] *eff[1]*eff[3]*eff[4] \
                   + P_hwp*apt_eff*eff[1]*emiss[1]  *eff[3]*eff[4]\
                   + P_aperture*(1.-apt_eff)*emiss[2]  *eff[3]*eff[4]\
                   + P_mirror*2.*eff[3]*emiss[3]  *eff[4]\
                   + P_1Kfilter*eff[4]*emiss[4] )

# Microstrip related efficiency and lenslet AR
# CMB, aperture
# HWP, emissivity, apt_eff, AR
# Aperture, 1-apt_eff
        P_load_nominal = 0.5*eff[6]*eff[5]  \
            * (P_cmb*apt_eff*eff[0]*emiss[0]   *eff[1]*eff[3]*eff[4]\
                   + P_hwp_nominal*apt_eff*eff[1]*emiss[1] *eff[3]*eff[4] \
                   + P_aperture_nominal*(1.-apt_eff)*emiss[2]   *eff[3]*eff[4]\
                   + P_mirror_nominal*2.*eff[3]*emiss[3]  *eff[4] \
                   + P_1Kfilter_nominal*eff[4]*emiss[4] )

        output[17,j]=P_load
        output[39,j]=P_cmb*apt_eff*eff[0]*emiss[0]*0.5*eff[6]*eff[5]

        poisson2 = np.sum( 2.*F_sum * g.h * freq * (freq[num-1]-freq[0])/float(num))
#        bunching2 = np.sum(F2_sum * (freq[num-1]-freq[0])/float(num))
#        bunching2 = np.sum(2.*F2_sum * (freq[num-1]-freq[0])/float(num))
        bunching2 = np.sum(2.*F_sum**2. * (freq[num-1]-freq[0])/float(num))
     
        NEP_ph_w = np.sqrt(poisson2 + bunching2)
        NEP_ph_wo = np.sqrt(poisson2)
        output[18,j]=NEP_ph_wo
        output[19,j]=NEP_ph_w
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; THERMAL NOISE

        if (P_load > Pmax): P_load_max = P_load_nominal
        if (P_load <= Pmax): 
            P_load_max = Pmax
            print 'maximum loading is taken for Jupiter', freq_c*1e-9, '%2.3f %2.3f %2.3f' % (P_load_nominal*1e12, P_load_max*1e12, Pmax*1e12)
        P_load_max = P_load_nominal
#;  print, '#############', T_bath
#        if (T_bath_nominal < 0.2): n_tc=1.   #  ; metal
#        if (T_bath_nominal > 0.2): n_tc=3.   #  ; phonon
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
        output[20,j]=NEP_th
        output[21,j]=T_bath
        output[22,j]=n_tc
        output[23,j]=T_bolo
        output[24,j]=gamma
        output[25,j]=X_fac
#        output[38,j]=X_fac*P_load_max/(T_bolo-T_bath)
        output[38,j]=Gave
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; READOUT NOISE
        R_bolo = 1.
        Vbias = np.sqrt((X_fac-1.)*R_bolo*P_load_max)
        NEP_readout = g.current_noise*Vbias
        output[37,j]=NEP_readout

#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; NEP_total     
        NEP_to_w = np.sqrt(NEP_ph_w**2. + NEP_th**2. + NEP_readout**2.)    #  ;*np.sqrt(2.)
        NEP_to_wo = np.sqrt(NEP_ph_wo**2. + NEP_th**2. + NEP_readout**2.)  #  ;*np.sqrt(2.)
        output[26,j]=NEP_to_w
        output[27,j]=NEP_to_wo
     
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
        NEQ_total_w = (np.sqrt(2.)*NET_totalper_w)/np.sqrt(float(num_pix))
        NEQ_total_wo = (np.sqrt(2.)*NET_totalper_wo)/np.sqrt(float(num_pix))
     
        output[28,j]=NET_photon_w
        output[29,j]=NET_photon_wo
        output[30,j]=NET_detect
        output[31,j]=NET_readout
        output[32,j]=NET_totalper_w
        output[33,j]=NET_totalper_wo
        output[34,j]=NEQ_total_w
        output[35,j]=NEQ_total_wo
        output[36,j]=NEP2NET_CMB(1.,freq,A_omega,eff_e2e)*1.e6
#  endfor

 
    filename = dir_out+'mappingspeed_'  \
        +str(int(freq_GHz))  \
        +'GHz_F'  \
        +'_T'+str(T_mir)+'K'
    print filename
    np.save(filename, output)

    return output
#end

#;####################################################################################################################################


def mapping_speed_dfix(d_pixel_mm, freq_GHz, Pmax, bandwidth, Npix, aperture_diameter_mm, emiss, eff, \
                       T_bath, T_mir, T_ape, T_hwp, T_1K, \
                       T_bath_nominal, T_mir_nominal, T_ape_nominal, T_hwp_nominal, T_1K_nominal, halfangle_edge_degs):

    num_elements = len(eff)
    eff_e2e_eae = 1.
    print 'system efficiency except aperture eff'
    for i in range(0, num_elements):
        eff_e2e_eae *= eff[i]
    print 'e2e except aperture eff',  eff_e2e_eae 
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
    sigma_0, Tedge_dB, apt_eff = aperture(D_lens,freq_c,halfangle_edge_degs)
    eff_e2e = eff_e2e_eae*apt_eff
    num_bolo = 2.*Npix
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; OPTICAL LOADING AND NEP_optical
 
#; define the frequency range
    freq_i = freq_c*(1.-bandwidth/2.)
    freq_f = freq_c*(1.+bandwidth/2.)
    freq = freq_i + np.arange(num)/float(num)*(freq_f-freq_i)
    A_omega = (g.c/freq)**2.
    output[0]=freq_c
    output[1]=freq_i
    output[2]=freq_f
    output[3]=bandwidth
    output[4]=np.mean(A_omega)
    output[5]=apt_eff

    output[6] = 0.
    output[7] = T_mir
    output[8] = Pmax
    output[9] = Npix
    output[10] = D_lens

#; integrate over the band
#;  CMB only
    BB = Brightness_BB(freq, g.Tcmb)
    F_cmb = flux(freq,g.Tcmb)
    P_cmb = np.sum( A_omega*BB * (freq[num-1]-freq[0])/float(num)) # ;*eff_opt*eff_spill*T_trans
    output[11]=P_cmb
    
#; HWP ony
    BB = Brightness_BB(freq, T_hwp)
    F_hwp = flux(freq,T_hwp)
    P_hwp = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*(1.-eff_spill)*T_trans
    output[12]=P_hwp

#; Nominal HWP temperature
    BB = Brightness_BB(freq, T_hwp_nominal)
    F_hwp_nominal = flux(freq,T_hwp_nominal)
    P_hwp_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #*eff_opt*(1.-eff_spill)*T_trans
      
#; Aperture only
    BB = Brightness_BB(freq, T_ape)
    F_aperture = flux(freq,T_ape)
    P_aperture = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*(1.-eff_spill)*T_trans
    output[12]=P_aperture

#; Nominal Aperture temperature
    BB = Brightness_BB(freq, T_ape_nominal)
    F_aperture_nominal = flux(freq,T_ape_nominal)
    P_aperture_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #*eff_opt*(1.-eff_spill)*T_trans
      
#; Mirror only
    BB = Brightness_BB(freq, T_mir)
    F_mirror = flux(freq,T_mir)
    P_mirror = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) # ;*eff_opt*emiss_mirror*2.*T_trans
    output[13]=P_mirror

#; Nominal Mirror temperature 
    BB = Brightness_BB(freq, T_mir_nominal)
    F_mirror_nominal = flux(freq,T_mir_nominal)
    P_mirror_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans

#; 1K filter only
    BB = Brightness_BB(freq, T_fil)
    F_1Kfilter = flux(freq,T_fil)
    P_1Kfilter = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans
    output[14]=P_1Kfilter

#; nominal 1K filter temperature
    BB = Brightness_BB(freq, T_1K_nominal)
    F_1Kfilter_nominal = flux(freq,T_1K_nominal)
    P_1Kfilter_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt*emiss_mirror*2.*T_trans
 
#; Pixel lens only
    BB = Brightness_BB(freq, T_bath)
    F_100mK = flux(freq, T_bath)
    P_100mK =  np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt
    output[15]=P_100mK

#; nominal Pixel lens temperature
    BB = Brightness_BB(freq, T_bath_nominal)
    F_100mK_nominal = flux(freq, T_bath_nominal)
    P_100mK_nominal = np.sum( A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt
 
#; Jupiter
    BB = Brightness_BB(freq, 1.86)                          #; 1.86K at 143GHz
    F_Jupiter = flux(freq,1.86)
    P_jupiter = 0.5*np.sum(A_omega*BB*(freq[num-1]-freq[0])/float(num)) #;*eff_opt
    output[16]=P_jupiter
 
    F_sum = eff[6]*eff[5] \
        *(F_cmb*apt_eff*eff[0]*emiss[0] *eff[1]*eff[3]*eff[4] \
          +F_hwp*apt_eff*eff[1]*emiss[1] *eff[3]*eff[4] \
          +F_aperture*(1.-apt_eff)*emiss[2] *eff[3]*eff[4] \
          +F_mirror*2.*eff[3]*emiss[3]  *emiss[4]\
          +F_1Kfilter*eff[4]*emiss[4] )

    F2_sum = (eff[6]*eff[5])**2  \
        * ( (F_cmb*apt_eff*eff[0]*emiss[0]  *eff[1]*eff[3]*eff[4])**2  \
                +(F_hwp*apt_eff*eff[1]*emiss[1] *eff[3]*eff[4])**2  \
                +(F_aperture*(1.-apt_eff)*emiss[2]   *eff[3]*eff[4])**2  \
                +(F_mirror*2.*eff[3]*emiss[3]  *eff[4])**2  \
                +(F_1Kfilter*eff[4]*emiss[4])**2 )              
    
    P_load = 0.5*eff[6]*eff[5]  \
        * (P_cmb*apt_eff*eff[0]*emiss[0] *eff[1]*eff[3]*eff[4] \
               + P_hwp*apt_eff*eff[1]*emiss[1]  *eff[3]*eff[4]\
               + P_aperture*(1.-apt_eff)*emiss[2]  *eff[3]*eff[4]\
               + P_mirror*2.*eff[3]*emiss[3]  *eff[4]\
               + P_1Kfilter*eff[4]*emiss[4] )

    P_load_nominal = 0.5*eff[6]*eff[5]  \
        * (P_cmb*apt_eff*eff[0]*emiss[0]   *eff[1]*eff[3]*eff[4]\
               + P_hwp_nominal*apt_eff*eff[1]*emiss[1] *eff[3]*eff[4] \
               + P_aperture_nominal*(1.-apt_eff)*emiss[2]   *eff[3]*eff[4]\
               + P_mirror_nominal*2.*eff[3]*emiss[3]  *eff[4] \
               + P_1Kfilter_nominal*eff[4]*emiss[4] )

    output[17]=P_load
    output[39]=P_cmb*apt_eff*eff[0]*emiss[0]*0.5*eff[6]*eff[5]

    poisson2 = np.sum( 2.*F_sum * g.h * freq * (freq[num-1]-freq[0])/float(num))
#        bunching2 = np.sum(F2_sum * (freq[num-1]-freq[0])/float(num))
#        bunching2 = np.sum(2.*F2_sum * (freq[num-1]-freq[0])/float(num))
    bunching2 = np.sum(2.*F_sum**2. * (freq[num-1]-freq[0])/float(num))
 
    NEP_ph_w = np.sqrt(poisson2 + bunching2)
    NEP_ph_wo = np.sqrt(poisson2)
    output[18]=NEP_ph_wo
    output[19]=NEP_ph_w
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; THERMAL NOISE

    if (P_load > Pmax): P_load_max = P_load_nominal
    if (P_load <= Pmax): 
        P_load_max = Pmax
        print 'maximum loading is taken for Jupiter', freq_c*1e-9, '%2.3f %2.3f %2.3f' % (P_load_nominal*1e12, P_load_max*1e12, Pmax*1e12)
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
    output[20]=NEP_th
    output[21]=T_bath
    output[22]=n_tc
    output[23]=T_bolo
    output[24]=gamma
    output[25]=X_fac
#        output[38,j]=X_fac*P_load_max/(T_bolo-T_bath)
    output[38]=Gave
#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; READOUT NOISE
    R_bolo = 1.
    Vbias = np.sqrt((X_fac-1.)*R_bolo*P_load_max)
    NEP_readout = g.current_noise*Vbias
    output[37]=NEP_readout

#;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#; NEP_total     
    NEP_to_w = np.sqrt(NEP_ph_w**2. + NEP_th**2. + NEP_readout**2.)    #  ;*np.sqrt(2.)
    NEP_to_wo = np.sqrt(NEP_ph_wo**2. + NEP_th**2. + NEP_readout**2.)  #  ;*np.sqrt(2.)
    output[26]=NEP_to_w
    output[27]=NEP_to_wo
 
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
 
    output[28]=NET_photon_w
    output[29]=NET_photon_wo
    output[30]=NET_detect
    output[31]=NET_readout
    output[32]=NET_totalper_w
    output[33]=NET_totalper_wo
    output[34]=NEQ_total_w
    output[35]=NEQ_total_wo
    output[36]=NEP2NET_CMB(1.,freq,A_omega,eff_e2e)*1.e6
#  endfor

    return output
#end

#;####################################################################################################################################


def mapping_speed_NdivNET2(dir_out):

    print '######################################################################'
    print ' aperture diameter: ', g.aperture_diameter_mm
    print ' mirror', g.T_mir_arr
    print ' aperture', g.T_ape
    print ' hwp', g.T_hwp
    print ' 1K', g.T_1K
    print ' bath temperature:', g.T_bath
    print 'halfangle_edge_degs', g.halfangle_edge_degs
    print ''
    print 'design temperatures are'
    print ' mirror', g.T_mir_nominal
    print ' aperture', g.T_ape_nominal
    print ' hwp', g.T_hwp_nominal
    print ' 1K', g.T_1K_nominal
    print ' bath', g.T_bath_nominal
    print '######################################################################'
    # This script runs for the NET for various bath temperature or
    # aperture temperature with the fixed bias voltage.
    
    # define basic parameter
    beam_ref = g.beam_ref             # str, i think this is the aparent aperture of the jupiter
    T_p = g.T_p                      # jupiter temperature
    ref_exp_beam = g.ref_exp_beam         # 1 degree at 100GHz
    
    freq_GHz_arr = g.freq_GHz_arr
    num_Tmir = len(g.T_mir_arr)
    dir_out = dir_out+'/'+str(int(g.T_bath*1e3))+'mK_Dapt'+str(int(g.aperture_diameter_mm))+'mm/data/'
    os.system('mkdir -p '+dir_out)
    
    print dir_out
    bandwidth = g.bandwidth
    wafer_size=g.wafer_size              # 8cm
    
    # CMB, HWP, Aperture, 2 Mirrors, 1K filter, Lens, Microstrip
    emiss = g.emiss                
    filename = dir_out+'emiss'
    np.save( filename+'.save', emiss)

    # CMB, 
    # HWP(AR), 
    # Aperture (not include in here), 
    # 2 mirrors, 1K filter 0.95 (goal value, Ade filter is 0.9-0.95), 
    # Lens AR, microstrip and antenna eff from Toki's PB value
    eff = g.eff   
    filename = dir_out+'eff'
    np.save(filename+'.save', eff)

    num_band = len(freq_GHz_arr)
    for i in range(0, num_band):
        for j in range(0,num_Tmir): 
            freq_GHz = freq_GHz_arr[i]
            theta = ref_exp_beam/180.*pi*(100./freq_GHz) / 2.
            beam_exp = 2.*pi*(1.-np.cos(theta))
            Pmax = Planet_dilution_T2P( freq_GHz, T_p, beam_ref, beam_exp, bandwidth[i])
            AOmega = (g.c/(freq_GHz*1e9))**2
            Tmax = abs(Brightness_dPdT( freq_GHz, g.Tcmb, AOmega, 0.68, bandwidth[i])**(-1) * Pmax)
            print ''
            print freq_GHz, ' GHz'
            print 'Jupiter beam incl. ', theta*radeg*2 ,' deg beam dilution', Pmax*1e12, ' pW' 
            print 'Jupiter beam in Kcmb', Tmax
      
            T_mir = g.T_mir_arr[j]
            mapping_speed( freq_GHz, Pmax, dir_out, \
                bandwidth[i], wafer_size, g.aperture_diameter_mm, \
                emiss, eff, \
                g.T_bath, T_mir, g.T_ape, g.T_hwp, g.T_1K, \
                g.T_bath_nominal, g.T_mir_nominal, g.T_ape_nominal, g.T_hwp_nominal, g.T_1K_nominal, g.halfangle_edge_degs)
