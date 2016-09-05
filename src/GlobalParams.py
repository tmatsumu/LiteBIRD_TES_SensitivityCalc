#####################################################################
#macbookpro
#dir1='/Users/tomotakematsumura/work/develop/POLARBEAR/Sensitivity/Cl_Nl/data_expectedCl/'
#bmode00
#dir1='/raid/users/tmatsumu/work/develop/PBI/Sensitivity/Cl_Nl/data_expectedCl/'
#dir1='/raid/users/tmatsumu/data_sim/Cls/wmap7BAOH0_params/'
dir_cl='/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20150208_ClTemplate/data/'
filename_lens = 'Class_planck2015_r0_cl_lensed.dat'
filename_prim = 'Class_planck2015_r1_nolens_cl.dat'
title = ''

fg_params = {}
fg_params['A_s'] = 0.012
fg_params['beta_s'] = -3.0
fg_params['m_s'] = -0.6
fg_params['nu_s'] = 90.e9

#fg_params['A_d'] = 0.004
fg_params['A_d'] = 0.1
fg_params['beta_d'] = 1.65
fg_params['m_d'] = -0.5
fg_params['nu_d'] = 90.e9

fg_params['ell_in0'] = 10.

nu1 = 60.e9
nu2 = 78.e9
nu3 = 100.e9
nu4 = 140.e9
nu5 = 195.e9
nu6 = 280.e9

r_input = 0.002

FWHM1 = 54.1
FWHM2 = 55.5
FWHM3 = 56.8
FWHM4 = 40.5
FWHM5 = 38.4
FWHM6 = 37.7

fsky = .5
noise_par = {}
noise_par['uKarcmin1'] = 15.72
noise_par['uKarcmin2'] = 9.86
noise_par['uKarcmin3'] = 7.06
noise_par['uKarcmin4'] = 5.59
noise_par['uKarcmin5'] = 4.70
noise_par['uKarcmin6'] = 5.69
noise_par['FWHM1'] = FWHM1
noise_par['FWHM2'] = FWHM2
noise_par['FWHM3'] = FWHM3
noise_par['FWHM4'] = FWHM4
noise_par['FWHM5'] = FWHM5
noise_par['FWHM6'] = FWHM6
noise_par['fsky1'] = fsky
noise_par['fsky2'] = fsky
noise_par['fsky3'] = fsky
noise_par['fsky4'] = fsky
noise_par['fsky5'] = fsky
noise_par['fsky6'] = fsky


fit_par = {}
fit_par = fg_params
fit_par['r'] = r_input

#dir_clfiles='/raid/users/tmatsumu/data_sim/Cls/wmap7BAOH0_params/'
#dir_out = '/home/tmatsumu/data_sim/LiteBIRD/Foreground/Analytic_template_subtraction/png/'
