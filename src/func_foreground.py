import numpy as np
import pylab as py
import lib_m as lib_m
import lib_foreground as lib_f
import lib_Clmanip as libcl
import GlobalParams as g
import sys

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   basic input parameters
pi = np.pi

#nu1 = 100.e9
#nu2 = 240.e9
#r_in = 0.1
#uKarcmin1 = 0.
#uKarcmin2 = 0.
#FWHM1 = 0.
#FWHM2 = 0.
#lmax = 300
#num_iter = 100
#option_unit = 'antenna'
##option_unit = 'thermo'
#option_fg = 'dust'

def func_foreground_2band2component(r_in,
								nu1, nu2,
								uKarcmin1, uKarcmin2,
								FWHM1, FWHM2,
								lmax, num_iter,
								option_fg, option_unit):
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   readin CMB Cl fits file, the cmb map is in thermodynamic unit by default
	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_L) != len(ell_P): print 'Something is wrong!'; sys.exit()
	ell_in = ell_L

	Dl_r1 = BBin_P
	Dl_lensing = BBin_L

	if option_unit == 'thermo':
	    A1 = 1.
	    A2 = 1.
	if option_unit == 'antenna':
	    A1 = lib_f.thermo2antenna_toDl(1.,nu1)
	    A2 = lib_f.thermo2antenna_toDl(1.,nu2)

# converting thermo unit to antenna unit if option_unit is antenna
	Dl_r1_nu1 = BBin_P*A1
	Dl_r1_nu2 = BBin_P*A2
	Dl_lensing_nu1 = BBin_L*A1
	Dl_lensing_nu2 = BBin_L*A2
	uKarcmin1 = uKarcmin1*A1
	uKarcmin2 = uKarcmin2*A2

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#   generalte the foreground model (no option means thermodynamic unit)

	ell_in, Cl_s1, Cl_d1 = lib_f.gen_Cl_Creminelli(ell_in,nu1,option_unit=option_unit)
	ell_in, Cl_s2, Cl_d2 = lib_f.gen_Cl_Creminelli(ell_in,nu2,option_unit=option_unit)

	if option_fg == 'synch': FG_synch = 1.; FG_dust = 0.
	if option_fg == 'dust': FG_synch = 0.; FG_dust = 1.
	if option_fg == 'both': FG_synch = 1.; FG_dust = 1.
	if option_fg == 'noFG': FG_synch = 0.; FG_dust = 0.
	Dl_s1 = lib_f.Cl2Dl(ell_in,Cl_s1*FG_synch)
	Dl_s2 = lib_f.Cl2Dl(ell_in,Cl_s2*FG_synch)
	Dl_d1 = lib_f.Cl2Dl(ell_in,Cl_d1*FG_dust)
	Dl_d2 = lib_f.Cl2Dl(ell_in,Cl_d2*FG_dust)

	# Cl inputs are assumed to thermodynamic unit. For antenna temperature, change to option_unit='antenna'
	ell_out, alpha_mean, alpha_std, r_mean, r_std, alpha_in = lib_f.Band2Comp2_FG(ell_in, r_in, \
																Dl_r1_nu1, Dl_r1_nu2, \
																Dl_lensing_nu1, Dl_lensing_nu2, \
										        				Dl_s1, Dl_s2, \
																Dl_d1, Dl_d2, \
																nu1, nu2, \
																uKarcmin1, uKarcmin2, \
																FWHM1, FWHM2, \
																lmax, num_iter, \
																option_unit)

	r_wm = lib_m.weight_mean_std(r_mean,r_std)
	alpha_wm = lib_m.weight_mean_std(alpha_mean,alpha_std)
	alphain_wm = np.mean(alpha_in)

	return r_wm, alpha_wm, alphain_wm


def func_foreground_2band2component_2dustcomp_v1(r_in,
								nu1, nu2,
								uKarcmin1, uKarcmin2,
								FWHM1, FWHM2,
								lmax, num_iter,
								option_fg, option_unit):
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   readin CMB Cl fits file, the cmb map is in thermodynamic unit by default
	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_L) != len(ell_P): print 'Something is wrong!'; sys.exit()
	ell_in = ell_L

	Dl_r1 = BBin_P
	Dl_lensing = BBin_L

	if option_unit == 'thermo':
	    A1 = 1.
	    A2 = 1.
	if option_unit == 'antenna':
	    A1 = lib_f.thermo2antenna_toDl(1.,nu1)
	    A2 = lib_f.thermo2antenna_toDl(1.,nu2)

# converting thermo unit to antenna unit if option_unit is antenna
	Dl_r1_nu1 = BBin_P*A1
	Dl_r1_nu2 = BBin_P*A2
	Dl_lensing_nu1 = BBin_L*A1
	Dl_lensing_nu2 = BBin_L*A2
	uKarcmin1 = uKarcmin1*A1
	uKarcmin2 = uKarcmin2*A2

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#   generalte the foreground model (no option means thermodynamic unit)

#	ell_in, Cl_s1, Cl_d1 = lib_f.gen_Cl_Creminelli(ell_in,nu1,option_unit=option_unit)
#	ell_in, Cl_s2, Cl_d2 = lib_f.gen_Cl_Creminelli(ell_in,nu2,option_unit=option_unit)
	ell_in, Cl_s1, Cl_d1, Cl_d11, Cl_d12, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu1,option_unit=option_unit)
	ell_in, Cl_s2, Cl_d2, Cl_d21, Cl_d22, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu2,option_unit=option_unit)

	if option_fg == 'synch': FG_synch = 1.; FG_dust = 0.
	if option_fg == 'dust': FG_synch = 0.; FG_dust = 1.
	if option_fg == 'both': FG_synch = 1.; FG_dust = 1.
	if option_fg == 'noFG': FG_synch = 0.; FG_dust = 0.
	Dl_s1 = lib_f.Cl2Dl(ell_in,Cl_s1*FG_synch)
	Dl_s2 = lib_f.Cl2Dl(ell_in,Cl_s2*FG_synch)
	Dl_d1 = lib_f.Cl2Dl(ell_in,Cl_d1*FG_dust)
	Dl_d2 = lib_f.Cl2Dl(ell_in,Cl_d2*FG_dust)

	# Cl inputs are assumed to thermodynamic unit. For antenna temperature, change to option_unit='antenna'
	ell_out, alpha_mean, alpha_std, r_mean, r_std, alpha_in = lib_f.Band2Comp2_FG(ell_in, r_in, \
																Dl_r1_nu1, Dl_r1_nu2, \
																Dl_lensing_nu1, Dl_lensing_nu2, \
										        				Dl_s1, Dl_s2, \
																Dl_d1, Dl_d2, \
																nu1, nu2, \
																uKarcmin1, uKarcmin2, \
																FWHM1, FWHM2, \
																lmax, num_iter, \
																option_unit)

	r_wm = lib_m.weight_mean_std(r_mean,r_std)
	alpha_wm = lib_m.weight_mean_std(alpha_mean,alpha_std)
	alphain_wm = np.mean(alpha_in)

	return r_wm, alpha_wm, alphain_wm

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   readin CMB Cl fits file, the cmb map is in thermodynamic unit by default
def func_foreground_3band2component(r_in,
								nu1, nu2, nu3,
								uKarcmin1, uKarcmin2, uKarcmin3,
								FWHM1, FWHM2, FWHM3,
								lmax, num_iter,
								option_unit):
	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_L) != len(ell_P): print 'Something is wrong!'; sys.exit()
	ell_in = ell_L

	Dl_r1 = BBin_P
	Dl_lensing = BBin_L

	if option_unit == 'thermo':
	    A1 = 1.
	    A2 = 1.
	    A3 = 1.
	if option_unit == 'antenna':
	    A1 = lib_f.thermo2antenna_toDl(1.,nu1)
	    A2 = lib_f.thermo2antenna_toDl(1.,nu2)
	    A3 = lib_f.thermo2antenna_toDl(1.,nu3)

# converting thermo unit to antenna unit if option_unit is antenna
	Dl_r1_nu1 = BBin_P*A1
	Dl_r1_nu2 = BBin_P*A2
	Dl_r1_nu3 = BBin_P*A3
	Dl_lensing_nu1 = BBin_L*A1
	Dl_lensing_nu2 = BBin_L*A2
	Dl_lensing_nu3 = BBin_L*A3
	uKarcmin1 = uKarcmin1*A1
	uKarcmin2 = uKarcmin2*A2
	uKarcmin3 = uKarcmin2*A3

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#   generalte the foreground model (no option means thermodynamic unit)

	ell_in, Cl_s1, Cl_d1 = lib_f.gen_Cl_Creminelli(ell_in,nu1,option_unit=option_unit)
	ell_in, Cl_s2, Cl_d2 = lib_f.gen_Cl_Creminelli(ell_in,nu2,option_unit=option_unit)
	ell_in, Cl_s3, Cl_d3 = lib_f.gen_Cl_Creminelli(ell_in,nu3,option_unit=option_unit)

	Dl_s1 = lib_f.Cl2Dl(ell_in,Cl_s1)
	Dl_s2 = lib_f.Cl2Dl(ell_in,Cl_s2)
	Dl_s3 = lib_f.Cl2Dl(ell_in,Cl_s3)
	Dl_d1 = lib_f.Cl2Dl(ell_in,Cl_d1)
	Dl_d2 = lib_f.Cl2Dl(ell_in,Cl_d2)
	Dl_d3 = lib_f.Cl2Dl(ell_in,Cl_d3)

	# Cl inputs are assumed to thermodynamic unit. For antenna temperature, change to option_unit='antenna'
	ell_out, alpha1_mean, alpha1_std, alpha3_mean, alpha3_std, r_mean, r_std, alpha_in \
		= lib_f.Band3Comp2_FG(ell_in, r_in, \
								Dl_r1_nu1, Dl_r1_nu2, Dl_r1_nu3, \
								Dl_lensing_nu1, Dl_lensing_nu2, Dl_lensing_nu3, \
								Dl_s1, Dl_s2, Dl_s3, \
								Dl_d1, Dl_d2, Dl_d3, \
								nu1, nu2, nu3, \
								uKarcmin1, uKarcmin2, uKarcmin3, \
								FWHM1, FWHM2, FWHM3, \
								lmax, num_iter, \
								option_unit)

	r_wm = lib_m.weight_mean_std(r_mean,r_std)
	alpha1_wm = lib_m.weight_mean_std(alpha1_mean,alpha1_std)
	alpha3_wm = lib_m.weight_mean_std(alpha3_mean,alpha3_std)

	alpha1in_wm = np.mean(alpha_in[0])
	alpha3in_wm = np.mean(alpha_in[1])

	return r_wm, alpha1_wm, alpha3_wm, alpha1in_wm, alpha3in_wm
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   readin CMB Cl fits file, the cmb map is in thermodynamic unit by default
def func_foreground_3band2component_2dustcomp_v1(r_in,
								nu1, nu2, nu3,
								uKarcmin1, uKarcmin2, uKarcmin3,
								FWHM1, FWHM2, FWHM3,
								lmax, num_iter,
								option_unit):
	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_L) != len(ell_P): print 'Something is wrong!'; sys.exit()
	ell_in = ell_L

	Dl_r1 = BBin_P
	Dl_lensing = BBin_L

	if option_unit == 'thermo':
	    A1 = 1.
	    A2 = 1.
	    A3 = 1.
	if option_unit == 'antenna':
	    A1 = lib_f.thermo2antenna_toDl(1.,nu1)
	    A2 = lib_f.thermo2antenna_toDl(1.,nu2)
	    A3 = lib_f.thermo2antenna_toDl(1.,nu3)

# converting thermo unit to antenna unit if option_unit is antenna
	Dl_r1_nu1 = BBin_P*A1
	Dl_r1_nu2 = BBin_P*A2
	Dl_r1_nu3 = BBin_P*A3
	Dl_lensing_nu1 = BBin_L*A1
	Dl_lensing_nu2 = BBin_L*A2
	Dl_lensing_nu3 = BBin_L*A3
	uKarcmin1 = uKarcmin1*A1
	uKarcmin2 = uKarcmin2*A2
	uKarcmin3 = uKarcmin2*A3

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#   generalte the foreground model (no option means thermodynamic unit)
	ell_in, Cl_s1, Cl_d1, Cl_d11, Cl_d12, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu1,option_unit=option_unit)
	ell_in, Cl_s2, Cl_d2, Cl_d21, Cl_d22, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu2,option_unit=option_unit)
	ell_in, Cl_s3, Cl_d3, Cl_d31, Cl_d32, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu3,option_unit=option_unit)
#	ell_in, Cl_s1, Cl_d1 = lib_f.gen_Cl_Creminelli(ell_in,nu1,option_unit=option_unit)
#	ell_in, Cl_s2, Cl_d2 = lib_f.gen_Cl_Creminelli(ell_in,nu2,option_unit=option_unit)
#	ell_in, Cl_s3, Cl_d3 = lib_f.gen_Cl_Creminelli(ell_in,nu3,option_unit=option_unit)

	Dl_s1 = lib_f.Cl2Dl(ell_in,Cl_s1)
	Dl_s2 = lib_f.Cl2Dl(ell_in,Cl_s2)
	Dl_s3 = lib_f.Cl2Dl(ell_in,Cl_s3)
	Dl_d1 = lib_f.Cl2Dl(ell_in,Cl_d1)
	Dl_d2 = lib_f.Cl2Dl(ell_in,Cl_d2)
	Dl_d3 = lib_f.Cl2Dl(ell_in,Cl_d3)

	# Cl inputs are assumed to thermodynamic unit. For antenna temperature, change to option_unit='antenna'
	ell_out, alpha1_mean, alpha1_std, alpha3_mean, alpha3_std, r_mean, r_std, alpha_in \
		= lib_f.Band3Comp2_FG(ell_in, r_in, \
								Dl_r1_nu1, Dl_r1_nu2, Dl_r1_nu3, \
								Dl_lensing_nu1, Dl_lensing_nu2, Dl_lensing_nu3, \
								Dl_s1, Dl_s2, Dl_s3, \
								Dl_d1, Dl_d2, Dl_d3, \
								nu1, nu2, nu3, \
								uKarcmin1, uKarcmin2, uKarcmin3, \
								FWHM1, FWHM2, FWHM3, \
								lmax, num_iter, \
								option_unit)

	r_wm = lib_m.weight_mean_std(r_mean,r_std)
	alpha1_wm = lib_m.weight_mean_std(alpha1_mean,alpha1_std)
	alpha3_wm = lib_m.weight_mean_std(alpha3_mean,alpha3_std)

	alpha1in_wm = np.mean(alpha_in[0])
	alpha3in_wm = np.mean(alpha_in[1])

	return r_wm, alpha1_wm, alpha3_wm, alpha1in_wm, alpha3in_wm
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   readin CMB Cl fits file, the cmb map is in thermodynamic unit by default
def func_noFG_3band(r_in,
					nu1, nu2, nu3,
					uKarcmin1, uKarcmin2, uKarcmin3,
					FWHM1, FWHM2, FWHM3,
					lmax, num_iter,
					option_unit):
	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_L) != len(ell_P): print 'Something is wrong!'; sys.exit()
	ell_in = ell_L

	Dl_r1 = BBin_P
	Dl_lensing = BBin_L

	if option_unit == 'thermo':
	    A1 = 1.
	    A2 = 1.
	    A3 = 1.
	if option_unit == 'antenna':
	    A1 = lib_f.thermo2antenna_toDl(1.,nu1)
	    A2 = lib_f.thermo2antenna_toDl(1.,nu2)
	    A3 = lib_f.thermo2antenna_toDl(1.,nu3)

# converting thermo unit to antenna unit if option_unit is antenna
	Dl_r1_nu1 = BBin_P*A1
	Dl_r1_nu2 = BBin_P*A2
	Dl_r1_nu3 = BBin_P*A3
	Dl_lensing_nu1 = BBin_L*A1
	Dl_lensing_nu2 = BBin_L*A2
	Dl_lensing_nu3 = BBin_L*A3
	uKarcmin1 = uKarcmin1*A1
	uKarcmin2 = uKarcmin2*A2
	uKarcmin3 = uKarcmin2*A3

	# Cl inputs are assumed to thermodynamic unit. For antenna temperature, change to option_unit='antenna'
	ell_out, alpha1_mean, alpha1_std, alpha3_mean, alpha3_std, r_mean, r_std \
		= lib_f.Band3Comp2_noFG(ell_in, r_in, \
								Dl_r1_nu1, Dl_r1_nu2, Dl_r1_nu3, \
								Dl_lensing_nu1, Dl_lensing_nu2, Dl_lensing_nu3, \
								nu1, nu2, nu3, \
								uKarcmin1, uKarcmin2, uKarcmin3, \
								FWHM1, FWHM2, FWHM3, \
								lmax, num_iter, \
								option_unit)

	r_wm = lib_m.weight_mean_std(r_mean,r_std)
	alpha1_wm = lib_m.weight_mean_std(alpha1_mean,alpha1_std)
	alpha3_wm = lib_m.weight_mean_std(alpha3_mean,alpha3_std)

	alpha1in_wm = np.mean(alpha_in[0])
	alpha3in_wm = np.mean(alpha_in[1])

	return r_wm, alpha1_wm, alpha3_wm, alpha1in_wm, alpha3in_wm
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   readin CMB Cl fits file, the cmb map is in thermodynamic unit by default
def func_foreground_6band2component(r_in,
								nu1, nu2, nu3, nu4, nu5, nu6,
								uKarcmin1, uKarcmin2, uKarcmin3, uKarcmin4, uKarcmin5, uKarcmin6,
								FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6,
								lmax, num_iter,
								option_unit):
	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_L) != len(ell_P): print 'Something is wrong!'; sys.exit()
	ell_in = ell_L

	Dl_r1 = BBin_P
	Dl_lensing = BBin_L

	if option_unit == 'thermo':
	    A1 = 1.
	    A2 = 1.
	    A3 = 1.
	    A4 = 1.
	    A5 = 1.
	    A6 = 1.
	if option_unit == 'antenna':
	    A1 = lib_f.thermo2antenna_toDl(1.,nu1)
	    A2 = lib_f.thermo2antenna_toDl(1.,nu2)
	    A3 = lib_f.thermo2antenna_toDl(1.,nu3)
	    A4 = lib_f.thermo2antenna_toDl(1.,nu4)
	    A5 = lib_f.thermo2antenna_toDl(1.,nu5)
	    A6 = lib_f.thermo2antenna_toDl(1.,nu6)

# converting thermo unit to antenna unit if option_unit is antenna
	Dl_r1_nu1 = BBin_P*A1
	Dl_r1_nu2 = BBin_P*A2
	Dl_r1_nu3 = BBin_P*A3
	Dl_r1_nu4 = BBin_P*A4
	Dl_r1_nu5 = BBin_P*A5
	Dl_r1_nu6 = BBin_P*A6
	Dl_lensing_nu1 = BBin_L*A1
	Dl_lensing_nu2 = BBin_L*A2
	Dl_lensing_nu3 = BBin_L*A3
	Dl_lensing_nu4 = BBin_L*A4
	Dl_lensing_nu5 = BBin_L*A5
	Dl_lensing_nu6 = BBin_L*A6
	uKarcmin1 = uKarcmin1*A1
	uKarcmin2 = uKarcmin2*A2
	uKarcmin3 = uKarcmin3*A3
	uKarcmin4 = uKarcmin4*A4
	uKarcmin5 = uKarcmin5*A5
	uKarcmin6 = uKarcmin6*A6

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#   generalte the foreground model (no option means thermodynamic unit)

	ell_in, Cl_s1, Cl_d1 = lib_f.gen_Cl_Creminelli(ell_in,nu1,option_unit=option_unit)
	ell_in, Cl_s2, Cl_d2 = lib_f.gen_Cl_Creminelli(ell_in,nu2,option_unit=option_unit)
	ell_in, Cl_s3, Cl_d3 = lib_f.gen_Cl_Creminelli(ell_in,nu3,option_unit=option_unit)
	ell_in, Cl_s4, Cl_d4 = lib_f.gen_Cl_Creminelli(ell_in,nu4,option_unit=option_unit)
	ell_in, Cl_s5, Cl_d5 = lib_f.gen_Cl_Creminelli(ell_in,nu5,option_unit=option_unit)
	ell_in, Cl_s6, Cl_d6 = lib_f.gen_Cl_Creminelli(ell_in,nu6,option_unit=option_unit)

	Dl_s1 = lib_f.Cl2Dl(ell_in,Cl_s1*1.)
	Dl_s2 = lib_f.Cl2Dl(ell_in,Cl_s2*1.)
	Dl_s3 = lib_f.Cl2Dl(ell_in,Cl_s3*1.)
	Dl_s4 = lib_f.Cl2Dl(ell_in,Cl_s4*1.)
	Dl_s5 = lib_f.Cl2Dl(ell_in,Cl_s5*1.)
	Dl_s6 = lib_f.Cl2Dl(ell_in,Cl_s6*1.)
	Dl_d1 = lib_f.Cl2Dl(ell_in,Cl_d1*1.)
	Dl_d2 = lib_f.Cl2Dl(ell_in,Cl_d2*1.)
	Dl_d3 = lib_f.Cl2Dl(ell_in,Cl_d3*1.)
	Dl_d4 = lib_f.Cl2Dl(ell_in,Cl_d4*1.)
	Dl_d5 = lib_f.Cl2Dl(ell_in,Cl_d5*1.)
	Dl_d6 = lib_f.Cl2Dl(ell_in,Cl_d6*1.)

	# Cl inputs are assumed to thermodynamic unit. For antenna temperature, change to option_unit='antenna'
	ell_out, alpha12_mean, alpha12_std, alpha65_mean, alpha65_std, r_mean, r_std, alpha_in \
		= lib_f.Band6Comp2_FG(ell_in, r_in, \
								Dl_r1_nu1, Dl_r1_nu2, Dl_r1_nu3, Dl_r1_nu4, Dl_r1_nu5, Dl_r1_nu6, \
								Dl_lensing_nu1, Dl_lensing_nu2, Dl_lensing_nu3, Dl_lensing_nu4, Dl_lensing_nu5, Dl_lensing_nu6, \
								Dl_s1, Dl_s2, Dl_s3, Dl_s4, Dl_s5, Dl_s6, \
								Dl_d1, Dl_d2, Dl_d3, Dl_d4, Dl_d5, Dl_d6, \
								nu1, nu2, nu3, nu4, nu5, nu6, \
								uKarcmin1, uKarcmin2, uKarcmin3, uKarcmin4, uKarcmin5, uKarcmin6, \
								FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6, \
								lmax, num_iter, \
								option_unit)

	r_wm = lib_m.weight_mean_std(r_mean,r_std)
	alpha12_wm = lib_m.weight_mean_std(alpha12_mean,alpha12_std)
	alpha65_wm = lib_m.weight_mean_std(alpha65_mean,alpha65_std)

	alpha12in_wm = np.mean(alpha_in[0])
	alpha65in_wm = np.mean(alpha_in[1])

	return r_wm, alpha12_wm, alpha65_wm, alpha12in_wm, alpha65in_wm

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   readin CMB Cl fits file, the cmb map is in thermodynamic unit by default
def func_foreground_6band2component_2dustcomp_v1(r_in,
								nu1, nu2, nu3, nu4, nu5, nu6,
								uKarcmin1, uKarcmin2, uKarcmin3, uKarcmin4, uKarcmin5, uKarcmin6,
								FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6,
								lmax, num_iter,
								option_unit, option_plot):
	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_L) != len(ell_P): print 'Something is wrong!'; sys.exit()
	ell_in = ell_L

	Dl_r1 = BBin_P
	Dl_lensing = BBin_L

	if option_unit == 'thermo':
	    A1 = 1.
	    A2 = 1.
	    A3 = 1.
	    A4 = 1.
	    A5 = 1.
	    A6 = 1.
	if option_unit == 'antenna':
	    A1 = lib_f.thermo2antenna(1.,nu1)
	    A2 = lib_f.thermo2antenna(1.,nu2)
	    A3 = lib_f.thermo2antenna(1.,nu3)
	    A4 = lib_f.thermo2antenna(1.,nu4)
	    A5 = lib_f.thermo2antenna(1.,nu5)
	    A6 = lib_f.thermo2antenna(1.,nu6)

# converting thermo unit to antenna unit if option_unit is antenna
	Dl_r1_nu1 = BBin_P*A1**2
	Dl_r1_nu2 = BBin_P*A2**2
	Dl_r1_nu3 = BBin_P*A3**2
	Dl_r1_nu4 = BBin_P*A4**2
	Dl_r1_nu5 = BBin_P*A5**2
	Dl_r1_nu6 = BBin_P*A6**2
	Dl_lensing_nu1 = BBin_L*A1**2
	Dl_lensing_nu2 = BBin_L*A2**2
	Dl_lensing_nu3 = BBin_L*A3**2
	Dl_lensing_nu4 = BBin_L*A4**2
	Dl_lensing_nu5 = BBin_L*A5**2
	Dl_lensing_nu6 = BBin_L*A6**2
	uKarcmin1 = uKarcmin1*A1
	uKarcmin2 = uKarcmin2*A2
	uKarcmin3 = uKarcmin3*A3
	uKarcmin4 = uKarcmin4*A4
	uKarcmin5 = uKarcmin5*A5
	uKarcmin6 = uKarcmin6*A6

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#   generalte the foreground model (no option means thermodynamic unit)

	ell_in, Cl_s1, Cl_d1, Cl_d11, Cl_d12, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu1,option_unit=option_unit)
	ell_in, Cl_s2, Cl_d2, Cl_d21, Cl_d22, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu2,option_unit=option_unit)
	ell_in, Cl_s3, Cl_d3, Cl_d31, Cl_d32, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu3,option_unit=option_unit)
	ell_in, Cl_s4, Cl_d4, Cl_d41, Cl_d42, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu4,option_unit=option_unit)
	ell_in, Cl_s5, Cl_d5, Cl_d51, Cl_d52, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu5,option_unit=option_unit)
	ell_in, Cl_s6, Cl_d6, Cl_d61, Cl_d62, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu6,option_unit=option_unit)


	if option_plot != '':
		num_nu = 100
		nu0_arr = np.logspace(10,np.log10(3000.e9),num_nu)
		Cl_rms2_PB = np.zeros(num_nu)
		Cl_rms2_LB = np.zeros(num_nu)
		Cl_rms2_Ds = np.zeros(num_nu)
		Cl_rms2_D1 = np.zeros(num_nu)
		Cl_rms2_D2 = np.zeros(num_nu)
		Cl_rms2_D0 = np.zeros(num_nu)
		Cl_rms2_S = np.zeros(num_nu)
		Cl_Pcmb = lib_f.Dl2Cl(ell_in,BBin_P)
		Cl_Lcmb = lib_f.Dl2Cl(ell_in,BBin_L)
		for i in range(num_nu):
			ell0, Cl_s0, Cl_d0, Cl_d01, Cl_d02, Cl_dsingle = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu0_arr[i],option_unit='thermo')
			sigma_b_in = 0.
			Cl_rms2_PB[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_Pcmb * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_LB[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_Lcmb * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_Ds[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_dsingle * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D1[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d01 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D2[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d02 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D0[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d0 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_S[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_s0 * np.exp(-ell_in**2*sigma_b_in**2))

		py.figure(figsize=(15,8))
		py.subplot(121)

		py.plot(nu0_arr*1.e-9,Cl_rms2_S,'r--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D0,'g--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D1,'b--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D2,'c--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_Ds,'m--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_S*1.e-4,'r')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D0*1.e-4,'g')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D1*1.e-4,'b')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D2*1.e-4,'c')
		py.plot(nu0_arr*1.e-9,Cl_rms2_Ds*1.e-4,'m')
		py.plot(nu0_arr*1.e-9, Cl_rms2_PB*r_in)
		py.plot(nu0_arr*1.e-9, Cl_rms2_LB)
		py.loglog()
		py.xticks( color = 'k', size = 15)
		py.yticks( color = 'k', size = 15)
		py.ylim([1e-7,1e4])
		py.xlim([np.min(nu0_arr*1e-9),np.max(nu0_arr*1e-9)])
		py.ylabel('$C_{rms}$ [$\mu K_{CMB}^{2}$]',fontsize=15)
#		py.ylabel('$C_{l=2}$ [$\mu$K$^{2}_{'+option_unit+'}$]',fontsize=15)
		py.xlabel('Frequency [GHz]',fontsize=15)

		num_nu = 100
		nu0_arr = np.logspace(10,np.log10(3000.e9),num_nu)
		Cl_rms2_PB = np.zeros(num_nu)
		Cl_rms2_LB = np.zeros(num_nu)
		Cl_rms2_Ds = np.zeros(num_nu)
		Cl_rms2_D1 = np.zeros(num_nu)
		Cl_rms2_D2 = np.zeros(num_nu)
		Cl_rms2_D0 = np.zeros(num_nu)
		Cl_rms2_S = np.zeros(num_nu)
		Cl_Pcmb = lib_f.Dl2Cl(ell_in,BBin_P)
		Cl_Lcmb = lib_f.Dl2Cl(ell_in,BBin_L)
		for i in range(num_nu):
			ell0, Cl_s0, Cl_d0, Cl_d01, Cl_d02, Cl_dsingle = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu0_arr[i],option_unit='antenna')
			sigma_b_in = 0.
			Cl_rms2_PB[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_Pcmb * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_LB[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_Lcmb * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_Ds[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_dsingle * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D1[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d01 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D2[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d02 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D0[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d0 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_S[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_s0 * np.exp(-ell_in**2*sigma_b_in**2))
		py.subplot(122)
		py.plot(nu0_arr*1.e-9,Cl_rms2_S,'r--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D0,'g--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D1,'b--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D2,'c--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_Ds,'m--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_S*1.e-4,'r')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D0*1.e-4,'g')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D1*1.e-4,'b')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D2*1.e-4,'c')
		py.plot(nu0_arr*1.e-9,Cl_rms2_Ds*1.e-4,'m')
		py.plot(nu0_arr*1.e-9, lib_f.thermo2antenna(Cl_rms2_PB,nu0_arr)*r_in)
		py.plot(nu0_arr*1.e-9, lib_f.thermo2antenna(Cl_rms2_LB,nu0_arr))
		py.loglog()
		py.xticks( color = 'k', size = 15)
		py.yticks( color = 'k', size = 15)
		py.ylim([1e-7,1e4])
		py.xlim([np.min(nu0_arr*1e-9),np.max(nu0_arr*1e-9)])
		py.ylabel('$C_{rms}$ [$\mu K_{antenna}^{2}$]',fontsize=15)
		py.xlabel('Frequency [GHz]',fontsize=15)
		print option_plot
		py.savefig(option_plot)
#		py.show()
		py.clf()
#		sys.exit()

	Dl_s1 = lib_f.Cl2Dl(ell_in,Cl_s1*1.)
	Dl_s2 = lib_f.Cl2Dl(ell_in,Cl_s2*1.)
	Dl_s3 = lib_f.Cl2Dl(ell_in,Cl_s3*1.)
	Dl_s4 = lib_f.Cl2Dl(ell_in,Cl_s4*1.)
	Dl_s5 = lib_f.Cl2Dl(ell_in,Cl_s5*1.)
	Dl_s6 = lib_f.Cl2Dl(ell_in,Cl_s6*1.)
	Dl_d1 = lib_f.Cl2Dl(ell_in,Cl_d1*1.)
	Dl_d2 = lib_f.Cl2Dl(ell_in,Cl_d2*1.)
	Dl_d3 = lib_f.Cl2Dl(ell_in,Cl_d3*1.)
	Dl_d4 = lib_f.Cl2Dl(ell_in,Cl_d4*1.)
	Dl_d5 = lib_f.Cl2Dl(ell_in,Cl_d5*1.)
	Dl_d6 = lib_f.Cl2Dl(ell_in,Cl_d6*1.)

	# Cl inputs are assumed to thermodynamic unit. For antenna temperature, change to option_unit='antenna'
	ell_out, alpha12_mean, alpha12_std, alpha65_mean, alpha65_std, r_mean, r_std, alpha_in \
		= lib_f.Band6Comp2_FG(ell_in, r_in, \
								Dl_r1_nu1, Dl_r1_nu2, Dl_r1_nu3, Dl_r1_nu4, Dl_r1_nu5, Dl_r1_nu6, \
								Dl_lensing_nu1, Dl_lensing_nu2, Dl_lensing_nu3, Dl_lensing_nu4, Dl_lensing_nu5, Dl_lensing_nu6, \
								Dl_s1, Dl_s2, Dl_s3, Dl_s4, Dl_s5, Dl_s6, \
								Dl_d1, Dl_d2, Dl_d3, Dl_d4, Dl_d5, Dl_d6, \
								nu1, nu2, nu3, nu4, nu5, nu6, \
								uKarcmin1, uKarcmin2, uKarcmin3, uKarcmin4, uKarcmin5, uKarcmin6, \
								FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6, \
								lmax, num_iter, \
								option_unit)

	r_wm = lib_m.weight_mean_std(r_mean,r_std)
	alpha12_wm = lib_m.weight_mean_std(alpha12_mean,alpha12_std)
	alpha65_wm = lib_m.weight_mean_std(alpha65_mean,alpha65_std)

	alpha12in_wm = np.mean(alpha_in[0])
	alpha65in_wm = np.mean(alpha_in[1])

	return r_wm, alpha12_wm, alpha65_wm, alpha12in_wm, alpha65in_wm
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   readin CMB Cl fits file, the cmb map is in thermodynamic unit by default
def func_foreground_6band2component_2dustcomp_v2(r_in,
								nu1, nu2, nu3, nu4, nu5, nu6,
								uKarcmin1, uKarcmin2, uKarcmin3, uKarcmin4, uKarcmin5, uKarcmin6,
								FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6,
								lmax, num_iter,
								option_unit, option_plot, factor_fg_s, factor_fg_d):
	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_L) != len(ell_P): print 'Something is wrong!'; sys.exit()
	ell_in = ell_L

	Dl_r1 = BBin_P
	Dl_lensing = BBin_L

	if option_unit == 'thermo':
	    A1 = 1.
	    A2 = 1.
	    A3 = 1.
	    A4 = 1.
	    A5 = 1.
	    A6 = 1.
	if option_unit == 'antenna':
	    A1 = lib_f.thermo2antenna(1.,nu1)
	    A2 = lib_f.thermo2antenna(1.,nu2)
	    A3 = lib_f.thermo2antenna(1.,nu3)
	    A4 = lib_f.thermo2antenna(1.,nu4)
	    A5 = lib_f.thermo2antenna(1.,nu5)
	    A6 = lib_f.thermo2antenna(1.,nu6)

# converting thermo unit to antenna unit if option_unit is antenna
	Dl_r1_nu1 = BBin_P*A1**2
	Dl_r1_nu2 = BBin_P*A2**2
	Dl_r1_nu3 = BBin_P*A3**2
	Dl_r1_nu4 = BBin_P*A4**2
	Dl_r1_nu5 = BBin_P*A5**2
	Dl_r1_nu6 = BBin_P*A6**2
	Dl_lensing_nu1 = BBin_L*A1**2
	Dl_lensing_nu2 = BBin_L*A2**2
	Dl_lensing_nu3 = BBin_L*A3**2
	Dl_lensing_nu4 = BBin_L*A4**2
	Dl_lensing_nu5 = BBin_L*A5**2
	Dl_lensing_nu6 = BBin_L*A6**2
	uKarcmin1 = uKarcmin1*A1
	uKarcmin2 = uKarcmin2*A2
	uKarcmin3 = uKarcmin3*A3
	uKarcmin4 = uKarcmin4*A4
	uKarcmin5 = uKarcmin5*A5
	uKarcmin6 = uKarcmin6*A6

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#   generalte the foreground model (no option means thermodynamic unit)

	ell_in, Cl_s1, Cl_d1, Cl_d11, Cl_d12, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu1,option_unit=option_unit)
	ell_in, Cl_s2, Cl_d2, Cl_d21, Cl_d22, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu2,option_unit=option_unit)
	ell_in, Cl_s3, Cl_d3, Cl_d31, Cl_d32, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu3,option_unit=option_unit)
	ell_in, Cl_s4, Cl_d4, Cl_d41, Cl_d42, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu4,option_unit=option_unit)
	ell_in, Cl_s5, Cl_d5, Cl_d51, Cl_d52, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu5,option_unit=option_unit)
	ell_in, Cl_s6, Cl_d6, Cl_d61, Cl_d62, Cl_tmp = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu6,option_unit=option_unit)


#	nu0_arr = np.logspace(10,np.log10(3000.e9),100)
#	ell0, Cl_s0, Cl_d0, Cl_d01, Cl_d02, Cl_dsingle = lib_f.gen_Cl_Creminellibased_2dustcomp(2,nu0_arr,option_unit='thermo')
#
#	if option_plot != '':
#		py.figure(figsize=(15,8))
#		py.subplot(121)
#		py.plot(nu0_arr*1.e-9,Cl_s0,'r--')
#		py.plot(nu0_arr*1.e-9,Cl_d0,'g--')
#		py.plot(nu0_arr*1.e-9,Cl_d01,'b--')
#		py.plot(nu0_arr*1.e-9,Cl_d02,'c--')
#		py.plot(nu0_arr*1.e-9,Cl_dsingle,'m--')
#		py.plot(nu0_arr*1.e-9,Cl_s0*1.e-4,'r')
#		py.plot(nu0_arr*1.e-9,Cl_d0*1.e-4,'g')
#		py.plot(nu0_arr*1.e-9,Cl_d01*1.e-4,'b')
#		py.plot(nu0_arr*1.e-9,Cl_d02*1.e-4,'c')
#		py.plot(nu0_arr*1.e-9,Cl_dsingle*1.e-4,'m')
#		py.plot(nu0_arr*1.e-9, lib_f.Dl2Cl(ell_in[0],BBin_P[0]*r_in*np.ones(len(nu0_arr))))
#		py.plot(nu0_arr*1.e-9, lib_f.Dl2Cl(ell_in[0],BBin_L[0]*np.ones(len(nu0_arr))))
#		py.loglog()
#		py.xticks( color = 'k', size = 15)
#		py.yticks( color = 'k', size = 15)
#		py.ylim([1e-7,1e4])
#		py.xlim([np.min(nu0_arr*1e-9),np.max(nu0_arr*1e-9)])
#		py.ylabel('$C_{l=2}$ [$\mu$K$^2_{'+option_unit+'}$]',fontsize=15)
#		py.xlabel('Frequency [GHz]',fontsize=15)
#
#		ell0, Cl_s0_ant, Cl_d0_ant, Cl_d01_ant, Cl_d02_ant, Cl_dsingle_ant = lib_f.gen_Cl_Creminellibased_2dustcomp(2,nu0_arr,option_unit='antenna')
#		py.subplot(122)
#		py.plot(nu0_arr*1.e-9,Cl_s0_ant,'r--')
#		py.plot(nu0_arr*1.e-9,Cl_d0_ant,'g--')
#		py.plot(nu0_arr*1.e-9,Cl_d01_ant,'b--')
#		py.plot(nu0_arr*1.e-9,Cl_d02_ant,'c--')
#		py.plot(nu0_arr*1.e-9,Cl_dsingle_ant,'m--')
#		py.plot(nu0_arr*1.e-9,Cl_s0_ant*1.e-4,'r')
#		py.plot(nu0_arr*1.e-9,Cl_d0_ant*1.e-4,'g')
##		py.plot(nu0_arr*1.e-9,Cl_d02_ant*1.e-4,'c')
#		py.plot(nu0_arr*1.e-9,Cl_dsingle_ant*1.e-4,'m')
#		py.plot(nu0_arr*1.e-9, lib_f.Dl2Cl(ell_in[0],BBin_P[0]*r_in*lib_f.thermo2antenna_toDl(1.,nu0_arr)))
#		py.plot(nu0_arr*1.e-9, lib_f.Dl2Cl(ell_in[0],BBin_L[0]*lib_f.thermo2antenna_toDl(1.,nu0_arr)))
#		py.loglog()
#		py.xticks( color = 'k', size = 15)
##		py.ylim([1e-7,1e4])
#		py.xlim([np.min(nu0_arr*1e-9),np.max(nu0_arr*1e-9)])
#		py.ylabel('$C_{l=2}$ [$\mu$K$^2_{antenna}$]',fontsize=15)
#		py.xlabel('Frequency [GHz]',fontsize=15)
#		py.savefig(option_plot)
##		py.show()
#		py.clf()
#		sys.exit()
	if option_plot != '':
		num_nu = 100
		nu0_arr = np.logspace(10,np.log10(3000.e9),num_nu)
		Cl_rms2_PB = np.zeros(num_nu)
		Cl_rms2_LB = np.zeros(num_nu)
		Cl_rms2_Ds = np.zeros(num_nu)
		Cl_rms2_D1 = np.zeros(num_nu)
		Cl_rms2_D2 = np.zeros(num_nu)
		Cl_rms2_D0 = np.zeros(num_nu)
		Cl_rms2_S = np.zeros(num_nu)
		Cl_Pcmb = lib_f.Dl2Cl(ell_in,BBin_P)
		Cl_Lcmb = lib_f.Dl2Cl(ell_in,BBin_L)
		for i in range(num_nu):
			ell0, Cl_s0, Cl_d0, Cl_d01, Cl_d02, Cl_dsingle = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu0_arr[i],option_unit='thermo')
			sigma_b_in = 0.
			Cl_rms2_PB[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_Pcmb * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_LB[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_Lcmb * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_Ds[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_dsingle * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D1[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d01 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D2[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d02 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D0[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d0 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_S[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_s0 * np.exp(-ell_in**2*sigma_b_in**2))

		py.figure(figsize=(15,8))
		py.subplot(121)

		py.plot(nu0_arr*1.e-9,Cl_rms2_S,'r--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D0,'g--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D1,'b--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D2,'c--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_Ds,'m--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_S*1.e-4,'r')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D0*1.e-4,'g')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D1*1.e-4,'b')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D2*1.e-4,'c')
		py.plot(nu0_arr*1.e-9,Cl_rms2_Ds*1.e-4,'m')
		py.plot(nu0_arr*1.e-9, Cl_rms2_PB*r_in)
		py.plot(nu0_arr*1.e-9, Cl_rms2_LB)
		py.loglog()
		py.xticks( color = 'k', size = 15)
		py.yticks( color = 'k', size = 15)
		py.ylim([1e-7,1e4])
		py.xlim([np.min(nu0_arr*1e-9),np.max(nu0_arr*1e-9)])
		py.ylabel('$C_{rms}$ [$\mu K_{CMB}^{2}$]',fontsize=15)
#		py.ylabel('$C_{l=2}$ [$\mu$K$^{2}_{'+option_unit+'}$]',fontsize=15)
		py.xlabel('Frequency [GHz]',fontsize=15)

		num_nu = 100
		nu0_arr = np.logspace(10,np.log10(3000.e9),num_nu)
		Cl_rms2_PB = np.zeros(num_nu)
		Cl_rms2_LB = np.zeros(num_nu)
		Cl_rms2_Ds = np.zeros(num_nu)
		Cl_rms2_D1 = np.zeros(num_nu)
		Cl_rms2_D2 = np.zeros(num_nu)
		Cl_rms2_D0 = np.zeros(num_nu)
		Cl_rms2_S = np.zeros(num_nu)
		Cl_Pcmb = lib_f.Dl2Cl(ell_in,BBin_P)
		Cl_Lcmb = lib_f.Dl2Cl(ell_in,BBin_L)
		for i in range(num_nu):
			ell0, Cl_s0, Cl_d0, Cl_d01, Cl_d02, Cl_dsingle = lib_f.gen_Cl_Creminellibased_2dustcomp(ell_in,nu0_arr[i],option_unit='antenna')
			sigma_b_in = 0.
			Cl_rms2_PB[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_Pcmb * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_LB[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_Lcmb * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_Ds[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_dsingle * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D1[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d01 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D2[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d02 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_D0[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_d0 * np.exp(-ell_in**2*sigma_b_in**2))
			Cl_rms2_S[i] = (1./(4.*pi))*np.sum( (2.*ell_in+1.) *Cl_s0 * np.exp(-ell_in**2*sigma_b_in**2))
		py.subplot(122)
		py.plot(nu0_arr*1.e-9,Cl_rms2_S,'r--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D0,'g--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D1,'b--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D2,'c--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_Ds,'m--')
		py.plot(nu0_arr*1.e-9,Cl_rms2_S*1.e-4,'r')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D0*1.e-4,'g')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D1*1.e-4,'b')
		py.plot(nu0_arr*1.e-9,Cl_rms2_D2*1.e-4,'c')
		py.plot(nu0_arr*1.e-9,Cl_rms2_Ds*1.e-4,'m')
		py.plot(nu0_arr*1.e-9, Cl_rms2_PB*r_in)
		py.plot(nu0_arr*1.e-9, Cl_rms2_LB)
		py.loglog()
		py.xticks( color = 'k', size = 15)
		py.yticks( color = 'k', size = 15)
		py.ylim([1e-7,1e4])
		py.xlim([np.min(nu0_arr*1e-9),np.max(nu0_arr*1e-9)])
		py.ylabel('$C_{rms}$ [$\mu K_{antenna}^{2}$]',fontsize=15)
		py.xlabel('Frequency [GHz]',fontsize=15)
		print option_plot
		py.savefig(option_plot)
#		py.show()
		py.clf()

	Dl_s1 = lib_f.Cl2Dl(ell_in,Cl_s1*factor_fg_s)
	Dl_s2 = lib_f.Cl2Dl(ell_in,Cl_s2*factor_fg_s)
	Dl_s3 = lib_f.Cl2Dl(ell_in,Cl_s3*factor_fg_s)
	Dl_s4 = lib_f.Cl2Dl(ell_in,Cl_s4*factor_fg_s)
	Dl_s5 = lib_f.Cl2Dl(ell_in,Cl_s5*factor_fg_s)
	Dl_s6 = lib_f.Cl2Dl(ell_in,Cl_s6*factor_fg_s)
	Dl_d1 = lib_f.Cl2Dl(ell_in,Cl_d1*factor_fg_d)
	Dl_d2 = lib_f.Cl2Dl(ell_in,Cl_d2*factor_fg_d)
	Dl_d3 = lib_f.Cl2Dl(ell_in,Cl_d3*factor_fg_d)
	Dl_d4 = lib_f.Cl2Dl(ell_in,Cl_d4*factor_fg_d)
	Dl_d5 = lib_f.Cl2Dl(ell_in,Cl_d5*factor_fg_d)
	Dl_d6 = lib_f.Cl2Dl(ell_in,Cl_d6*factor_fg_d)

	# Cl inputs are assumed to thermodynamic unit. For antenna temperature, change to option_unit='antenna'
	ell_out, alpha12_mean, alpha12_std, alpha65_mean, alpha65_std, r_mean, r_std, alpha_in \
		= lib_f.Band6Comp2_FG(ell_in, r_in, \
								Dl_r1_nu1, Dl_r1_nu2, Dl_r1_nu3, Dl_r1_nu4, Dl_r1_nu5, Dl_r1_nu6, \
								Dl_lensing_nu1, Dl_lensing_nu2, Dl_lensing_nu3, Dl_lensing_nu4, Dl_lensing_nu5, Dl_lensing_nu6, \
								Dl_s1, Dl_s2, Dl_s3, Dl_s4, Dl_s5, Dl_s6, \
								Dl_d1, Dl_d2, Dl_d3, Dl_d4, Dl_d5, Dl_d6, \
								nu1, nu2, nu3, nu4, nu5, nu6, \
								uKarcmin1, uKarcmin2, uKarcmin3, uKarcmin4, uKarcmin5, uKarcmin6, \
								FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6, \
								lmax, num_iter, \
								option_unit)

	r_wm = lib_m.weight_mean_std(r_mean,r_std)
	alpha12_wm = lib_m.weight_mean_std(alpha12_mean,alpha12_std)
	alpha65_wm = lib_m.weight_mean_std(alpha65_mean,alpha65_std)

	alpha12in_wm = np.mean(alpha_in[0])
	alpha65in_wm = np.mean(alpha_in[1])

	return r_wm, alpha12_wm, alpha65_wm, alpha12in_wm, alpha65in_wm
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   readin CMB Cl fits file, the cmb map is in thermodynamic unit by default
def func_foreground_6band2component_fgcorr(r_in,
								nu1, nu2, nu3, nu4, nu5, nu6,
								uKarcmin1, uKarcmin2, uKarcmin3, uKarcmin4, uKarcmin5, uKarcmin6,
								FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6,
								lmax, num_iter,
								option_unit):
	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_L) != len(ell_P): print 'Something is wrong!'; sys.exit()
	ell_in = ell_L

	Dl_r1 = BBin_P
	Dl_lensing = BBin_L

	if option_unit == 'thermo':
	    A1 = 1.
	    A2 = 1.
	    A3 = 1.
	    A4 = 1.
	    A5 = 1.
	    A6 = 1.
	if option_unit == 'antenna':
	    A1 = lib_f.thermo2antenna_toDl(1.,nu1)
	    A2 = lib_f.thermo2antenna_toDl(1.,nu2)
	    A3 = lib_f.thermo2antenna_toDl(1.,nu3)
	    A4 = lib_f.thermo2antenna_toDl(1.,nu4)
	    A5 = lib_f.thermo2antenna_toDl(1.,nu5)
	    A6 = lib_f.thermo2antenna_toDl(1.,nu6)

# converting thermo unit to antenna unit if option_unit is antenna
	Dl_r1_nu1 = BBin_P*A1
	Dl_r1_nu2 = BBin_P*A2
	Dl_r1_nu3 = BBin_P*A3
	Dl_r1_nu4 = BBin_P*A4
	Dl_r1_nu5 = BBin_P*A5
	Dl_r1_nu6 = BBin_P*A6
	Dl_lensing_nu1 = BBin_L*A1
	Dl_lensing_nu2 = BBin_L*A2
	Dl_lensing_nu3 = BBin_L*A3
	Dl_lensing_nu4 = BBin_L*A4
	Dl_lensing_nu5 = BBin_L*A5
	Dl_lensing_nu6 = BBin_L*A6
	uKarcmin1 = uKarcmin1*A1
	uKarcmin2 = uKarcmin2*A2
	uKarcmin3 = uKarcmin3*A3
	uKarcmin4 = uKarcmin4*A4
	uKarcmin5 = uKarcmin5*A5
	uKarcmin6 = uKarcmin6*A6

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#   generalte the foreground model (no option means thermodynamic unit)

	ell_in, Cl_s1, Cl_d1 = lib_f.gen_Cl_Creminelli(ell_in,nu1,option_unit=option_unit)
	ell_in, Cl_s2, Cl_d2 = lib_f.gen_Cl_Creminelli(ell_in,nu2,option_unit=option_unit)
	ell_in, Cl_s3, Cl_d3 = lib_f.gen_Cl_Creminelli(ell_in,nu3,option_unit=option_unit)
	ell_in, Cl_s4, Cl_d4 = lib_f.gen_Cl_Creminelli(ell_in,nu4,option_unit=option_unit)
	ell_in, Cl_s5, Cl_d5 = lib_f.gen_Cl_Creminelli(ell_in,nu5,option_unit=option_unit)
	ell_in, Cl_s6, Cl_d6 = lib_f.gen_Cl_Creminelli(ell_in,nu6,option_unit=option_unit)

	Dl_s1 = lib_f.Cl2Dl(ell_in,Cl_s1*1.)
	Dl_s2 = lib_f.Cl2Dl(ell_in,Cl_s2*1.)
	Dl_s3 = lib_f.Cl2Dl(ell_in,Cl_s3*1.)
	Dl_s4 = lib_f.Cl2Dl(ell_in,Cl_s4*1.)
	Dl_s5 = lib_f.Cl2Dl(ell_in,Cl_s5*1.)
	Dl_s6 = lib_f.Cl2Dl(ell_in,Cl_s6*1.)
	Dl_d1 = lib_f.Cl2Dl(ell_in,Cl_d1*1.)
	Dl_d2 = lib_f.Cl2Dl(ell_in,Cl_d2*1.)
	Dl_d3 = lib_f.Cl2Dl(ell_in,Cl_d3*1.)
	Dl_d4 = lib_f.Cl2Dl(ell_in,Cl_d4*1.)
	Dl_d5 = lib_f.Cl2Dl(ell_in,Cl_d5*1.)
	Dl_d6 = lib_f.Cl2Dl(ell_in,Cl_d6*1.)

	# Cl inputs are assumed to thermodynamic unit. For antenna temperature, change to option_unit='antenna'
	ell_out, alpha12_mean, alpha12_std, alpha65_mean, alpha65_std, r_mean, r_std, alpha_in \
		= lib_f.Band6Comp2_FGcorr(ell_in, r_in, \
								Dl_r1_nu1, Dl_r1_nu2, Dl_r1_nu3, Dl_r1_nu4, Dl_r1_nu5, Dl_r1_nu6, \
								Dl_lensing_nu1, Dl_lensing_nu2, Dl_lensing_nu3, Dl_lensing_nu4, Dl_lensing_nu5, Dl_lensing_nu6, \
								Dl_s1, Dl_s2, Dl_s3, Dl_s4, Dl_s5, Dl_s6, \
								Dl_d1, Dl_d2, Dl_d3, Dl_d4, Dl_d5, Dl_d6, \
								nu1, nu2, nu3, nu4, nu5, nu6, \
								uKarcmin1, uKarcmin2, uKarcmin3, uKarcmin4, uKarcmin5, uKarcmin6, \
								FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6, \
								lmax, num_iter, \
								option_unit)

	r_wm = lib_m.weight_mean_std(r_mean,r_std)
	alpha12_wm = lib_m.weight_mean_std(alpha12_mean,alpha12_std)
	alpha65_wm = lib_m.weight_mean_std(alpha65_mean,alpha65_std)

	alpha12in_wm = np.mean(alpha_in[0])
	alpha65in_wm = np.mean(alpha_in[1])

	return r_wm, alpha12_wm, alpha65_wm, alpha12in_wm, alpha65in_wm
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   readin CMB Cl fits file, the cmb map is in thermodynamic unit by default
def func_noFG_6band(r_in,
					nu1, nu2, nu3, nu4, nu5, nu6,
					uKarcmin1, uKarcmin2, uKarcmin3, uKarcmin4, uKarcmin5, uKarcmin6,
					FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6,
					lmax, num_iter,
					option_unit):
	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_L) != len(ell_P): print 'Something is wrong!'; sys.exit()
	ell_in = ell_L

	Dl_r1 = BBin_P
	Dl_lensing = BBin_L

	if option_unit == 'thermo':
	    A1 = 1.
	    A2 = 1.
	    A3 = 1.
	    A4 = 1.
	    A5 = 1.
	    A6 = 1.
	if option_unit == 'antenna':
	    A1 = lib_f.thermo2antenna_toDl(1.,nu1)
	    A2 = lib_f.thermo2antenna_toDl(1.,nu2)
	    A3 = lib_f.thermo2antenna_toDl(1.,nu3)
	    A4 = lib_f.thermo2antenna_toDl(1.,nu4)
	    A5 = lib_f.thermo2antenna_toDl(1.,nu5)
	    A6 = lib_f.thermo2antenna_toDl(1.,nu6)

# converting thermo unit to antenna unit if option_unit is antenna
	Dl_r1_nu1 = BBin_P*A1
	Dl_r1_nu2 = BBin_P*A2
	Dl_r1_nu3 = BBin_P*A3
	Dl_r1_nu4 = BBin_P*A4
	Dl_r1_nu5 = BBin_P*A5
	Dl_r1_nu6 = BBin_P*A6
	Dl_lensing_nu1 = BBin_L*A1
	Dl_lensing_nu2 = BBin_L*A2
	Dl_lensing_nu3 = BBin_L*A3
	Dl_lensing_nu4 = BBin_L*A4
	Dl_lensing_nu5 = BBin_L*A5
	Dl_lensing_nu6 = BBin_L*A6
	uKarcmin1 = uKarcmin1*A1
	uKarcmin2 = uKarcmin2*A2
	uKarcmin3 = uKarcmin3*A3
	uKarcmin4 = uKarcmin4*A4
	uKarcmin5 = uKarcmin5*A5
	uKarcmin6 = uKarcmin6*A6

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#   generalte the foreground model (no option means thermodynamic unit)

	# Cl inputs are assumed to thermodynamic unit. For antenna temperature, change to option_unit='antenna'
	ell_out, r_mean, r_std \
		= lib_f.Band6_noFG(ell_in, r_in, \
						Dl_r1_nu1, Dl_r1_nu2, Dl_r1_nu3, Dl_r1_nu4, Dl_r1_nu5, Dl_r1_nu6, \
						Dl_lensing_nu1, Dl_lensing_nu2, Dl_lensing_nu3, Dl_lensing_nu4, Dl_lensing_nu5, Dl_lensing_nu6, \
						nu1, nu2, nu3, nu4, nu5, nu6, \
						uKarcmin1, uKarcmin2, uKarcmin3, uKarcmin4, uKarcmin5, uKarcmin6, \
						FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6, \
						lmax, num_iter, \
						option_unit)

	r_wm = lib_m.weight_mean_std(r_mean,r_std)

	return r_wm
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def plotCl_2band(nu1,nu2,uKarcmin1,uKarcmin2,FWHM1,FWHM2,dir_out,filename_out):

	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_P) != len(ell_L): print 'something is wrong!'
	ell_in = ell_P

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	ell_in, Cl_s1, Cl_d1 = lib_f.gen_Cl_Creminelli(ell_in,nu1,option_unit='thermo')
	ell_in, Cl_s2, Cl_d2 = lib_f.gen_Cl_Creminelli(ell_in,nu2,option_unit='thermo')

	Dl_s1 = lib_f.Cl2Dl(ell_in,Cl_s1)
	Dl_s2 = lib_f.Cl2Dl(ell_in,Cl_s2)
	Dl_d1 = lib_f.Cl2Dl(ell_in,Cl_d1)
	Dl_d2 = lib_f.Cl2Dl(ell_in,Cl_d2)

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	Cl_zero = np.zeros(len(ell_in))

	gen_Nl = libcl.gen_Nl(ell_in)
	prefact = gen_Nl.cal_prefact()
	gen_Nl.C_l = Cl_zero
	gen_Nl.fsky = 1.
	gen_Nl.prefact_option(True)
	gen_Nl.modeloss_option(False)
	gen_Nl.uKarcmin = uKarcmin1
	gen_Nl.FWHM = FWHM1
	gen_Nl.sigma_b()
	Nl_1 = gen_Nl.gen_KnoxdNl('noCV')

	gen_Nl = libcl.gen_Nl(ell_in)
	prefact = gen_Nl.cal_prefact()
	gen_Nl.C_l = Cl_zero
	gen_Nl.fsky = 1.
	gen_Nl.prefact_option(True)
	gen_Nl.modeloss_option(False)
	gen_Nl.uKarcmin = uKarcmin2
	gen_Nl.FWHM = FWHM2
	gen_Nl.sigma_b()
	Nl_2 = gen_Nl.gen_KnoxdNl('noCV')
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	py.figure()
	py.plot(ell_P,BBin_P*0.1,'m',linewidth=2)
	py.plot(ell_P,BBin_P*0.01,'m',linewidth=2)
	py.plot(ell_P,BBin_P*0.001,'m',linewidth=2)
	py.plot(ell_L,BBin_L,'b',linewidth=2)
	py.plot(ell_P, BBin_P*0.1+BBin_L,'g',linewidth=2)
	py.plot(ell_P, BBin_P*0.01+BBin_L,'g',linewidth=2)
	py.plot(ell_P, BBin_P*0.001+BBin_L,'g',linewidth=2)

	py.plot(ell_in,Dl_s1,'c--',linewidth=2,label='Synch')
	py.plot(ell_in,Dl_s2,'c-.',linewidth=2)

	py.plot(ell_in,Dl_d1,'r--',linewidth=2,label='Dust')
	py.plot(ell_in,Dl_d2,'r-.',linewidth=2)

	py.plot(ell_in,Nl_1,'m--',linewidth=2,label='N_{l1}')
	py.plot(ell_in,Nl_2,'m-.',linewidth=2,label='N_{l2}')

	py.xlabel('$l$', fontsize=17)
	py.ylabel('$l(l+1)/2\pi C_l$ [$\mu$K$^2$]', fontsize=17)
	py.xticks( color = 'k', size = 17)
	py.yticks( color = 'k', size = 17)
	py.title('2-band, $\\nu_1$=%1.0f, $\\nu_2$=%1.0f' % (nu1*1e-9, nu2*1e-9))
	py.legend(loc='best')

	py.loglog()

	py.savefig(dir_out+'/png/'+filename_out+'_Cl2band')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plotCl_3band(nu1,nu2,nu3,uKarcmin1,uKarcmin2,uKarcmin3,FWHM1,FWHM2,FWHM3,dir_out,filename_out):

	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_P) != len(ell_L): print 'something is wrong!'
	ell_in = ell_P

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	ell_in, Cl_s1, Cl_d1 = lib_f.gen_Cl_Creminelli(ell_in,nu1,option_unit='thermo')
	ell_in, Cl_s2, Cl_d2 = lib_f.gen_Cl_Creminelli(ell_in,nu2,option_unit='thermo')
	ell_in, Cl_s3, Cl_d3 = lib_f.gen_Cl_Creminelli(ell_in,nu3,option_unit='thermo')

	Dl_s1 = lib_f.Cl2Dl(ell_in,Cl_s1)
	Dl_s2 = lib_f.Cl2Dl(ell_in,Cl_s2)
	Dl_s3 = lib_f.Cl2Dl(ell_in,Cl_s3)
	Dl_d1 = lib_f.Cl2Dl(ell_in,Cl_d1)
	Dl_d2 = lib_f.Cl2Dl(ell_in,Cl_d2)
	Dl_d3 = lib_f.Cl2Dl(ell_in,Cl_d3)

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	Cl_zero = np.zeros(len(ell_in))

	gen_Nl = libcl.gen_Nl(ell_in)
	prefact = gen_Nl.cal_prefact()
	gen_Nl.C_l = Cl_zero
	gen_Nl.fsky = 1.
	gen_Nl.prefact_option(True)
	gen_Nl.modeloss_option(False)
	gen_Nl.uKarcmin = uKarcmin1
	gen_Nl.FWHM = FWHM1
	gen_Nl.sigma_b()
	Nl_1 = gen_Nl.gen_KnoxdNl('noCV')

	gen_Nl = libcl.gen_Nl(ell_in)
	prefact = gen_Nl.cal_prefact()
	gen_Nl.C_l = Cl_zero
	gen_Nl.fsky = 1.
	gen_Nl.prefact_option(True)
	gen_Nl.modeloss_option(False)
	gen_Nl.uKarcmin = uKarcmin2
	gen_Nl.FWHM = FWHM2
	gen_Nl.sigma_b()
	Nl_2 = gen_Nl.gen_KnoxdNl('noCV')

	gen_Nl = libcl.gen_Nl(ell_in)
	prefact = gen_Nl.cal_prefact()
	gen_Nl.C_l = Cl_zero
	gen_Nl.fsky = 1.
	gen_Nl.prefact_option(True)
	gen_Nl.modeloss_option(False)
	gen_Nl.uKarcmin = uKarcmin3
	gen_Nl.FWHM = FWHM3
	gen_Nl.sigma_b()
	Nl_3 = gen_Nl.gen_KnoxdNl('noCV')
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	py.figure()
	py.plot(ell_P,BBin_P*0.1,'m',linewidth=2)
	py.plot(ell_P,BBin_P*0.01,'m',linewidth=2)
	py.plot(ell_P,BBin_P*0.001,'m',linewidth=2)
	py.plot(ell_L,BBin_L,'b',linewidth=2)
	py.plot(ell_P, BBin_P*0.1+BBin_L,'g',linewidth=2)
	py.plot(ell_P, BBin_P*0.01+BBin_L,'g',linewidth=2)
	py.plot(ell_P, BBin_P*0.001+BBin_L,'g',linewidth=2)

	py.plot(ell_in,Dl_s1,'c--',linewidth=2,label='Synch')
	py.plot(ell_in,Dl_s2,'c-.',linewidth=2)
	py.plot(ell_in,Dl_s3,'c:',linewidth=2)

	py.plot(ell_in,Dl_d1,'r--',linewidth=2,label='Dust')
	py.plot(ell_in,Dl_d2,'r-.',linewidth=2)
	py.plot(ell_in,Dl_d3,'r:',linewidth=2)

	py.plot(ell_in,Nl_1,'m--',linewidth=2,label='$N_{l1}$')
	py.plot(ell_in,Nl_2,'m-.',linewidth=2,label='$N_{l2}$')
	py.plot(ell_in,Nl_3,'m:',linewidth=2,label='$N_{l3}$')

	py.xlabel('$l$', fontsize=17)
	py.ylabel('$l(l+1)/2\pi C_l$ [$\mu$K$^2$]', fontsize=17)
	py.xticks( color = 'k', size = 17)
	py.yticks( color = 'k', size = 17)
	py.title('3-band, $\\nu_1$=%1.0f, $\\nu_2$=%1.0f, $\\nu_3$=%1.0f' % (nu1*1e-9, nu2*1e-9, nu3*1e-9))
	py.legend(loc='best')

	py.loglog()

	py.savefig(dir_out+'/png/'+filename_out+'_Cl3band')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plotCl_6band(nu1,nu2,nu3,nu4,nu5,nu6, \
				uKarcmin1,uKarcmin2,uKarcmin3,uKarcmin4,uKarcmin5,uKarcmin6, \
				FWHM1,FWHM2,FWHM3,FWHM4,FWHM5,FWHM6, \
				dir_out,filename_out):
	print ''
	print ' func_foreground.py: plotCl_6band'
	print ''

	read_obj = libcl.read_cambdata()
	Cl_prim  = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_prim)
	Cl_lens = read_obj.read_cl_classdata_in_cambformat(g.dir_cl+'/'+g.filename_lens)

	ell_P = np.array(np.int_(Cl_prim['ell']))
	EEin_P = Cl_prim['EE']
	BBin_P = Cl_prim['BB']

	ell_L = np.array(np.int_(Cl_lens['ell']))
	EEin_L = Cl_lens['EE']
	BBin_L = Cl_lens['BB']

	if len(ell_P) != len(ell_L): print 'something is wrong!'
	ell_in = ell_P

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	ell_in, Cl_s1, Cl_d1 = lib_f.gen_Cl_Creminelli(ell_in,nu1,option_unit='thermo')
	ell_in, Cl_s2, Cl_d2 = lib_f.gen_Cl_Creminelli(ell_in,nu2,option_unit='thermo')
	ell_in, Cl_s3, Cl_d3 = lib_f.gen_Cl_Creminelli(ell_in,nu3,option_unit='thermo')
	ell_in, Cl_s4, Cl_d4 = lib_f.gen_Cl_Creminelli(ell_in,nu4,option_unit='thermo')
	ell_in, Cl_s5, Cl_d5 = lib_f.gen_Cl_Creminelli(ell_in,nu5,option_unit='thermo')
	ell_in, Cl_s6, Cl_d6 = lib_f.gen_Cl_Creminelli(ell_in,nu6,option_unit='thermo')

	Dl_s1 = lib_f.Cl2Dl(ell_in,Cl_s1)
	Dl_s2 = lib_f.Cl2Dl(ell_in,Cl_s2)
	Dl_s3 = lib_f.Cl2Dl(ell_in,Cl_s3)
	Dl_s4 = lib_f.Cl2Dl(ell_in,Cl_s4)
	Dl_s5 = lib_f.Cl2Dl(ell_in,Cl_s5)
	Dl_s6 = lib_f.Cl2Dl(ell_in,Cl_s6)
	Dl_d1 = lib_f.Cl2Dl(ell_in,Cl_d1)
	Dl_d2 = lib_f.Cl2Dl(ell_in,Cl_d2)
	Dl_d3 = lib_f.Cl2Dl(ell_in,Cl_d3)
	Dl_d4 = lib_f.Cl2Dl(ell_in,Cl_d4)
	Dl_d5 = lib_f.Cl2Dl(ell_in,Cl_d5)
	Dl_d6 = lib_f.Cl2Dl(ell_in,Cl_d6)

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	Cl_zero = np.zeros(len(ell_in))

	gen_Nl = libcl.gen_Nl(ell_in)
	prefact = gen_Nl.cal_prefact()
	gen_Nl.C_l = Cl_zero
	gen_Nl.fsky = 1.
	gen_Nl.prefact_option(True)
	gen_Nl.modeloss_option(False)
	gen_Nl.uKarcmin = uKarcmin1
	gen_Nl.FWHM = FWHM1
	gen_Nl.sigma_b()
	Nl_1 = gen_Nl.gen_KnoxdNl('noCV')

	gen_Nl = libcl.gen_Nl(ell_in)
	prefact = gen_Nl.cal_prefact()
	gen_Nl.C_l = Cl_zero
	gen_Nl.fsky = 1.
	gen_Nl.prefact_option(True)
	gen_Nl.modeloss_option(False)
	gen_Nl.uKarcmin = uKarcmin2
	gen_Nl.FWHM = FWHM2
	gen_Nl.sigma_b()
	Nl_2 = gen_Nl.gen_KnoxdNl('noCV')

	gen_Nl = libcl.gen_Nl(ell_in)
	prefact = gen_Nl.cal_prefact()
	gen_Nl.C_l = Cl_zero
	gen_Nl.fsky = 1.
	gen_Nl.prefact_option(True)
	gen_Nl.modeloss_option(False)
	gen_Nl.uKarcmin = uKarcmin3
	gen_Nl.FWHM = FWHM3
	gen_Nl.sigma_b()
	Nl_3 = gen_Nl.gen_KnoxdNl('noCV')

	gen_Nl = libcl.gen_Nl(ell_in)
	prefact = gen_Nl.cal_prefact()
	gen_Nl.C_l = Cl_zero
	gen_Nl.fsky = 1.
	gen_Nl.prefact_option(True)
	gen_Nl.modeloss_option(False)
	gen_Nl.uKarcmin = uKarcmin4
	gen_Nl.FWHM = FWHM4
	gen_Nl.sigma_b()
	Nl_4 = gen_Nl.gen_KnoxdNl('noCV')

	gen_Nl = libcl.gen_Nl(ell_in)
	prefact = gen_Nl.cal_prefact()
	gen_Nl.C_l = Cl_zero
	gen_Nl.fsky = 1.
	gen_Nl.prefact_option(True)
	gen_Nl.modeloss_option(False)
	gen_Nl.uKarcmin = uKarcmin5
	gen_Nl.FWHM = FWHM5
	gen_Nl.sigma_b()
	Nl_5 = gen_Nl.gen_KnoxdNl('noCV')

	gen_Nl = libcl.gen_Nl(ell_in)
	prefact = gen_Nl.cal_prefact()
	gen_Nl.C_l = Cl_zero
	gen_Nl.fsky = 1.
	gen_Nl.prefact_option(True)
	gen_Nl.modeloss_option(False)
	gen_Nl.uKarcmin = uKarcmin6
	gen_Nl.FWHM = FWHM6
	gen_Nl.sigma_b()
	Nl_6 = gen_Nl.gen_KnoxdNl('noCV')
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	py.figure()
	py.plot(ell_P,BBin_P*0.1,'m',linewidth=2)
	py.plot(ell_P,BBin_P*0.01,'m',linewidth=2)
	py.plot(ell_P,BBin_P*0.001,'m',linewidth=2)
	py.plot(ell_L,BBin_L,'b',linewidth=2)
	py.plot(ell_P, BBin_P*0.1+BBin_L,'g',linewidth=2)
	py.plot(ell_P, BBin_P*0.01+BBin_L,'g',linewidth=2)
	py.plot(ell_P, BBin_P*0.001+BBin_L,'g',linewidth=2)

	py.plot(ell_in,Dl_s1,'c--',linewidth=2,label='Synch')
	py.plot(ell_in,Dl_s2,'c-.',linewidth=2)
	py.plot(ell_in,Dl_s3,'c:',linewidth=2)
	py.plot(ell_in,Dl_s4,'c--',linewidth=2)
	py.plot(ell_in,Dl_s5,'c-.',linewidth=2)
	py.plot(ell_in,Dl_s6,'c:',linewidth=2)

	py.plot(ell_in,Dl_d1,'r--',linewidth=2,label='Dust')
	py.plot(ell_in,Dl_d2,'r-.',linewidth=2)
	py.plot(ell_in,Dl_d3,'r:',linewidth=2)
	py.plot(ell_in,Dl_d4,'r-',linewidth=2)
	py.plot(ell_in,Dl_d5,'r-.',linewidth=2)
	py.plot(ell_in,Dl_d6,'r:',linewidth=2)

	py.plot(ell_in,Nl_1,'m--',linewidth=2,label='$N_{l1}$')
	py.plot(ell_in,Nl_2,'m-.',linewidth=2,label='$N_{l2}$')
	py.plot(ell_in,Nl_3,'m:',linewidth=2,label='$N_{l3}$')
	py.plot(ell_in,Nl_4,'m--',linewidth=2,label='$N_{l4}$')
	py.plot(ell_in,Nl_5,'m-.',linewidth=2,label='$N_{l5}$')
	py.plot(ell_in,Nl_6,'m:',linewidth=2,label='$N_{l6}$')

	py.xlabel('$l$', fontsize=17)
	py.ylabel('$l(l+1)/2\pi C_l$ [$\mu$K$^2$]', fontsize=17)
	py.xticks( color = 'k', size = 17)
	py.yticks( color = 'k', size = 17)
	py.title('6-band, $\\nu_1$=%1.0f, $\\nu_2$=%1.0f, $\\nu_3$=%1.0f, $\\nu_4$=%1.0f, $\\nu_5$=%1.0f, $\\nu_6$=%1.0f' % (nu1*1e-9, nu2*1e-9, nu3*1e-9, nu4*1e-9, nu5*1e-9, nu6*1e-9))
	py.legend(loc='best')

	py.loglog()

	py.savefig(dir_out+'/png/'+filename_out+'_Cl6band')