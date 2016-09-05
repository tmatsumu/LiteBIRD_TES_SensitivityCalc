import numpy as np
import pylab as py
import lib_mappingspeed as lib_ms
import func_foreground as func
import global_par as g
import os
import sys

pi = np.pi
dir_out = '/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20150701_LB_Sensitivity/data/susceptability/'
#++++++++++++++++++++++++++++++++++++++++++++++++++++++
r_in = 0.
fsky = 0.53

FWHM1=FWHM2=FWHM3=FWHM4=FWHM5=FWHM6=0.
lmax = 20
num_iter = 200
option_unit = 'thermo'
option_plot = ''
num_monte = 20

#++++++++++++++++++++++++++++++++++++++++++++++++++++++

nu_arr_GHz = g.freq_GHz_arr
nu_arr = np.array(g.freq_GHz_arr) * 1.e9
num_band = len(nu_arr)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++

t_obs = g.Tmis_sec
margin = g.uKarcmin_total_margin

#++++++++++++++++++++++++++++++++++++++++++++++++++++++
num_multichroic = 33

#option_trial_arr = ['T_bath', 'T_mir', 'T_ape', 'T_hwp', 'T_1K', \
#				'T_bath_nominal', 'T_mir_nominal', 'T_ape_nominal', 'T_hwp_nominal', 'T_1K_nominal', \
#				'emiss_mir', 'eff_hwp', 'eff_filter', 'eff_lenslet', 'eff_det', 'margin']
#option_trial_arr = ['T_1K','T_bath_nominal', 'T_mir_nominal', 'T_ape_nominal', 'T_hwp_nominal', 'T_1K_nominal', \
#				'emiss_mir', 'eff_hwp', 'eff_filter', 'eff_lenslet', 'eff_det']
#option_trial_arr = ['emiss_mir', 'eff_hwp', 'eff_filter', 'eff_lenslet', 'eff_det']
#option_trial_arr = ['eff_det']
#option_trial_arr = ['T_bath','T_bath_nominal']
#option_trial_arr = ['T_bath', 'T_ape', 'T_hwp', 'T_bath_nominal', 'T_ape_nominal', 'T_hwp_nominal', 'eff_hwp', 'eff_filter', 'eff_lenslet', 'eff_det', 'margin']
option_trial_arr = ['margin']

for option_trial in option_trial_arr:
	filename_out= dir_out + '/multichroic'+str(num_multichroic)+'_'+option_trial
	num_trials = 10
	r_arr = np.zeros((2,num_trials))
	uKarcmin_tot = []
	for j in range(num_trials):
#	for j in range(0,1,1):
		print ''
		print j, '/', num_trials, option_trial, num_multichroic

		if num_multichroic ==1:
			d_pixel_mm = np.array([20.35441739, 13.33067464,12.0710959,  5.60113676, 7.3461988, 12.61276229])
			Npix = np.array([45,103,102, 602,19,43])

		if num_multichroic ==2:
			d_pixel_mm = np.array([21.30891992, 21.30891992, 12.19068749, 12.19068749, 6.6938434, 6.6938434])
			Npix = np.array([99,99,170, 170,184,184])

		if num_multichroic ==3:
			d_pixel_mm = np.array([18.1727542,18.1727542,18.1727542, 12.4198804,12.4198804,12.4198804])
			Npix = np.array([179,179,179, 127,127,127])

		if num_multichroic ==6:
			d_pixel_mm = 22.84097924*np.ones(6)
			Npix = 151*np.ones(6)

		if num_multichroic ==33:
			d_pixel_mm = np.array([18.,18.,18., 12.,12.,12.])
			Npix = np.array([152,152,152, 185,185,185])

		if num_multichroic ==22:
			d_pixel_mm = np.array([21., 21., 12., 12., 6.7, 6.7])
			Npix = np.array([109,109, 188,188, 204,204])

		Pmax = 0.
		freq_GHz_arr = [60.,78.,100.,140.,195.,280.]
		bandwidth = [0.23, 0.23, 0.23, 0.3, 0.3, 0.3]
		aperture_diameter_mm = 400.

		# sky, HWP-AR, aperture, baffle?, filter, lenslet, detector
		emiss = [1.,0.1,1.,0.005,0.1,0.,1.]
		if option_trial == 'emiss_mir':
			emiss_mir_arr = np.hstack([0.005, np.linspace(0.001,0.05,num_trials-1)] )
			syspar_arr = emiss_mir_arr
			emiss = [1.,0.1,1.,emiss_mir_arr[j],0.1,0.,1.]
			print emiss_mir_arr[j]

#		eff = [1.,0.98,1.,1.,0.95,0.99,0.73]
		eff = [1.,0.98,1.,1.,0.95,0.9,0.809]
		if option_trial == 'eff_hwp':
			eff_hwp_arr = np.hstack([0.98, np.linspace(0.6,1.,num_trials-1)] )
			syspar_arr = eff_hwp_arr
			eff = [1.,eff_hwp_arr[j],1.,1.,0.95,0.99,0.73]
			print eff_hwp_arr[j]

		if option_trial == 'eff_filter':
			eff_filter_arr = np.hstack([0.95, np.linspace(0.6,1.,num_trials-1)])
			syspar_arr = eff_filter_arr
			eff = [1.,0.98,1.,1.,eff_filter_arr[j],0.99,0.73]
			print eff_filter_arr[j]

		if option_trial == 'eff_lenslet':
			eff_lenslet_arr = np.hstack([0.99, np.linspace(0.6,1.,num_trials-1)])
			syspar_arr = eff_lenslet_arr
			eff = [1.,0.98,1.,1.,0.95,eff_lenslet_arr[j],0.73]
			print eff_lenslet_arr[j]

		if option_trial == 'eff_det':
			eff_det_arr = np.hstack([0.73, np.linspace(0.4,0.9,num_trials-1)])
			syspar_arr = eff_det_arr
			eff = [1.,0.98,1.,1.,0.95,0.99,eff_det_arr[j]]
			print eff_det_arr[j]

		T_bath = 0.1
		if option_trial == 'T_bath':
			T_bath_arr = np.hstack([0.1, np.linspace(0.09,0.2,num_trials-1)])
			T_bath = T_bath_arr[j]
			syspar_arr = T_bath_arr
			print T_bath

		T_mir = 4.
		if option_trial == 'T_mir':
			T_mir_arr = np.hstack([4, np.linspace(2,100,num_trials-1)])
			T_mir = T_mir_arr[j]
			syspar_arr = T_mir_arr

		T_ape = 4.
		if option_trial == 'T_ape':
			T_ape_arr = np.hstack([4, np.linspace(2,10,num_trials-1)])
			T_ape = T_ape_arr[j]
			syspar_arr = T_ape_arr

		T_hwp = 4.
		if option_trial == 'T_hwp':
			T_hwp_arr = np.hstack([4, np.linspace(2,10,num_trials-1)])
			T_hwp = T_hwp_arr[j]
			syspar_arr = T_hwp_arr

		T_1K = 1.
		if option_trial == 'T_1K':
			T_1K_arr = np.hstack([1, np.linspace(0.5,10,num_trials-1)])
			T_1K = T_1K_arr[j]
			syspar_arr = T_1K_arr

		T_bath_nominal = 0.1
		if option_trial == 'T_bath_nominal':
			T_bath_nominal_arr = np.hstack([0.1, np.linspace(0.09,0.2,num_trials-1)])
			T_bath_nominal = T_bath_nominal_arr[j]
			syspar_arr = T_bath_nominal_arr
			T_1K = T_bath_nominal

		T_mir_nominal = 4.
		if option_trial == 'T_mir_nominal':
			T_mir_nominal_arr = np.hstack([4., np.linspace(2,100,num_trials-1)])
			T_mir_nominal = T_mir_nominal_arr[j]
			syspar_arr = T_mir_nominal_arr
			T_mir = T_mir_nominal

		T_ape_nominal = 4.
		if option_trial == 'T_ape_nominal':
			T_ape_nominal_arr = np.hstack([4., np.linspace(2,10,num_trials-1)])
			T_ape_nominal = T_ape_nominal_arr[j]
			syspar_arr = T_ape_nominal_arr
			T_ape = T_ape_nominal

		T_hwp_nominal = 4.
		if option_trial == 'T_hwp_nominal':
			T_hwp_nominal_arr = np.hstack([4., np.linspace(2,10,num_trials-1)])
			T_hwp_nominal = T_hwp_nominal_arr[j]
			syspar_arr = T_hwp_nominal_arr
			T_hwp = T_hwp_nominal

		T_1K_nominal = 1.
		if option_trial == 'T_1K_nominal':
			T_1K_nominal_arr = np.hstack([1, np.linspace(0.5,10,num_trials-1)])
			T_1K_nominal = T_1K_nominal_arr[j]
			syspar_arr = T_1K_nominal_arr
			T_1K = T_1K_nominal

		margin = g.uKarcmin_total_margin
		if option_trial == 'margin':
			margin_arr = np.hstack([g.uKarcmin_total_margin, np.linspace(0.2,1,num_trials-1)])
			margin = margin_arr[j]
			syspar_arr = margin_arr

		halfangle_edge_degs = 8.13
		#++++++++++++++++++++++++++++++++++++++++++++++++++++++

		uKarcmin_arr = np.zeros(num_band)
		NET_arr = []
		Npix_arr = []
		for i in range(num_band):
			out = lib_ms.mapping_speed_dfix(d_pixel_mm[i], freq_GHz_arr[i], Pmax, bandwidth[i], Npix[i], \
								aperture_diameter_mm, emiss, eff, \
		                       	T_bath, T_mir, T_ape, T_hwp, T_1K, \
		                       	T_bath_nominal, T_mir_nominal, T_ape_nominal, T_hwp_nominal, T_1K_nominal, \
		                       	halfangle_edge_degs)
			NET_arr.append(out[32])
			Npix_arr.append(out[9])

		#++++++++++++++++++++++++++++++++++++++++++++++++++++++

		uKarcmin_arr = 10800./pi * np.sqrt(8.*pi*np.array(NET_arr)**2/t_obs/(2.*np.array(Npix_arr)))/margin
		uKarcmin_tot.append( 1./np.sqrt(np.sum(1./uKarcmin_arr**2)) )
		#++++++++++++++++++++++++++++++++++++++++++++++++++++++
		print NET_arr
		print uKarcmin_arr
		r_wm_i_ = []
		for i in range(num_monte):
			r_wm_i, alpha12_wm_i, alpha65_wm_i, alpha12in_wm_i, alpha65in_wm_i = func.func_foreground_6band2component_2dustcomp_v1( r_in,
																				nu_arr[0], nu_arr[1], nu_arr[2], nu_arr[3], nu_arr[4], nu_arr[5], 
																				uKarcmin_arr[0], uKarcmin_arr[1], uKarcmin_arr[2], uKarcmin_arr[3], uKarcmin_arr[4], uKarcmin_arr[5], 
																				FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6,
																				lmax, num_iter,
																				option_unit,option_plot)
			r_wm_i_.append(r_wm_i[1]/np.sqrt(fsky))
		r_arr[0,j] = np.mean(np.array(r_wm_i_))
		r_arr[1,j] = np.std(np.array(r_wm_i_))
		print r_arr[1,j]
	#++++++++++++++++++++++++++++++++++++++++++++++++++++++
	np.savez(filename_out, delta_r=r_arr, syspar_arr=syspar_arr, uKarcmin_tot=uKarcmin_tot, 
		num_multichroic=num_multichroic, option_trial=option_trial, 
		help='delta_r, syspar_arr, uKarcmin_tot, option_trial, num_multichroic')
	print ''
	print syspar_arr, r_arr[1,:]
	py.errorbar(syspar_arr, r_arr[0,:], r_arr[1,:], fmt='o')
	py.title(option_trial)
#py.show()
