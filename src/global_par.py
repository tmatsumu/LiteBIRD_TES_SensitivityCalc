import numpy as np

# basic parameters
Tcmb = 2.725 
h = 6.626068e-34
k_b = 1.3806503e-23
c = 2.99792458e8 

# jupiter info
beam_ref = 2.481e-8  
T_p = 173.5
ref_exp_beam = 45./60.

r_in = 0.

# LB specificaitons

#option_proposal = 'JAXA_MDR'
#option_proposal = 'US_MO'
#option_proposal = 'US_MO_LFT'
#option_proposal = 'US_MO_HFT'
option_proposal = 'US_CSR_LFT'

if option_proposal == 'JAXA_MDR':
	freq_GHz_arr = [60.,78.,100.,140.,195.,280.]
	bandwidth_arr = [0.23, 0.23, 0.23, 0.3, 0.3, 0.3]
	wafer_num = [8,8,8,5,5,5]
	freq_break = 101.

if option_proposal == 'US_MO':
	freq_GHz_arr = [40., 50., 60., 68., 78., 89., 100., 119., 140., 166., 195., 235., 280., 337., 402.]
	bandwidth_arr = [   .3,  .3,  .23, .23, .23, .23,  .23,   .3,   .3,   .3,   .3,   .3,   .3,   .3,  .23]
	wafer_num = [8,8,8,8,8,8,5,5,5,5,5,5,1,1,1]
	freq_break = 99.

if option_proposal == 'JAXA_MDR_perwafer':
	freq_GHz_arr = [60.,78.,100.,140.,195.,280.]
	bandwidth_arr = [0.23, 0.23, 0.23, 0.3, 0.3, 0.3]
	wafer_num = [1,1,1,1,1,1]
	freq_break = 101.

if option_proposal == 'US_MO_LFT':
	freq_GHz_arr = [40., 50., 60., 68., 78., 89.,       100., 119., 140., 166., 195., 235.]
	bandwidth_arr = [   .3,  .3,  .23, .23, .23, .23,    .23,   .3,   .3,   .3,   .3,   .3]
	wafer_num = [4,4,4,4,4,4,   3,2,3,2,3,2]
	pixDmm = [18.,18.,18.,18.,18.,18.,   12.,12.,12.,12.,12.,12.]
	stupid_offset = [0.,0.,0.,0.,0.,0.,   42.,28.,42.,28.,42.,28.]
       	# sky, HWP, aperture, M1, M2, filter, lenslet, det
	emiss_arr = [1., .03, 1.0, 0.005, 0.005, 0.02, 1., 1.000]
        ref_arr = [0., 0.1,  0.,       0.001, 0.001, 0.1,    0.1,   0.3]
        
        #-------------- HWP parameters ------------------------
       	# assume a sapphire based HWP
        n_hwp = 3.4
	losstan_hwp = 5.e-5
	d_hwp = 27.e-3  # This number is the maximum number that comes from the sapphire 9 stack

        

	Fnum = 3.5
#	halfangle_edge_degs = np.arcsin(0.5/Fnum) *180./np.pi  # F/#=3.5
	halfangle_edge_degs = np.arctan(0.5/Fnum) *180./np.pi  # F/#=3.5
	beamwaistfactor = 2.6
	wafer_size = 80.e-3
	aperture_diameter_mm = 400.


	# we are not using a lens and so the parameters below aren't relevant
	n_lens = 1.
	losstan_lens = 5.e-5
	d_lens = 20.e-3

	# assume a PP (a base material for the metal mesh)
	n_filter = 1.5
	losstan_filter = 2.3e-4
	d_filter = 5.e-3  # just by seeing a Giampaolo MM filter

	# absorptive loss of the mirror
	mirror_abs = 0.005  # skin depth absorption from Charlie

if option_proposal == 'US_MO_HFT':
	freq_GHz_arr = [280., 337., 402.]
	bandwidth_arr = [ .3,   .3,  .23]
	wafer_num = [1,1,1]
	pixDmm = [5.2,4.5,3.8]
	stupid_offset = [-8.,-20.,52.]

	#           sky, HWP, aperture, L1, L2, filter, lenslet, det
	emiss_arr = [1., 0.03, 1.,    0.03, 0.03, 0.1,   0.,   0.]
	#         sky, HWP, aperture, L1, L2, filter, lenslet, det
#	eff_arr = [1., 0.9, 1.,      0.98, 0.98, 0.9, 1, 0.7]
	ref_arr = [0., 0.02, 0.,      0.02, 0.02, 0.02, 0., 0.3]
#	halfangle_edge_degs = 12.804   # F/#=2.2
	Fnum = 2.2
#	halfangle_edge_degs = np.arcsin(0.5/Fnum) *180./np.pi
	halfangle_edge_degs = np.arctan(0.5/Fnum) *180./np.pi
	beamwaistfactor = 3.1
	wafer_size = 30.e-3
	aperture_diameter_mm = 200.

	# assume a silicon based HWP
	n_hwp = 3.4
	losstan_hwp = 5.e-5
	d_hwp = 3.99e-3  # This number is the maximum number that comes from the sapphire 3 stack

	# assume a silicon lens
	n_lens = 3.4
	losstan_lens = 5.e-5
	d_lens = 20.e-3 # this is a thickness from the current optical design. This can be thinner but the concern is if this holds the launch vibration.

	# assume a PP (a base material for the metal mesh)
	n_filter = 1.5
	losstan_filter = 2.3e-4
	d_filter = 5.e-3  # just by seeing a Giampaolo MM filter

if option_proposal == 'US_CSR_LFT':
	freq_GHz_arr = [40., 50., 60., 68., 78., 89.,       100., 119., 140., 166., 195., 235.]
	bandwidth_arr = [   .3,  .3,  .23, .23, .23, .23,    .23,   .3,   .3,   .3,   .3,   .3]
	wafer_num = [4,4,4,4,4,4,   3,2,3,2,3,2]
	pixDmm = [18.,18.,18.,18.,18.,18.,   12.,12.,12.,12.,12.,12.]
        stupid_offset = [38.,38.,38.,38.,38.,38.,   -32.,-46.,-32.,-46.,-32.,-46.] 
	# sky, HWP, aperture, M1, M2, filter, lenslet, det
        T_ref_arr = [2.725, 5., 5., 5., 5., 0.1, 0.1, 0.1] # temperature that element reflects to
	emiss_arr = [1., .03, 1.0, 0.005, 0.005, 0.05, 1., 0.]
	ref_arr = [0., 0.08,  0., 0., 0., 0.05, 0.05, 0.32] 
        eff_arr = [1.,0.,0.,0.,0.,0.,0.,0.]
        NEP_tot_arr = [5.47, 6.1, 5.84, 6.13, 6.56, 6.93, 7., 8.47, 8.89, 9.14, 9.23, 9.06]
     
	Fnum = 3.5
	halfangle_edge_degs = np.arcsin(0.5/Fnum) *180./np.pi  # F/#=3.5
	beamwaistfactor = 2.6
	wafer_size = 80.e-3
	aperture_diameter_mm = 400.

	# assume a sapphire based HWP
	n_hwp = 3.24
	losstan_hwp = 5.e-5
	d_hwp = 27.e-3  # This number is the maximum number that comes from the sapphire 9 stack

	# we are not using a lens and so the parameters below aren't relevant
	n_lens = 3.4
	losstan_lens = 5.e-5
	d_lens = 9.e-3

	# assume a PP (a base material for the metal mesh)
	n_filter = 1.5
	losstan_filter = 2.3e-4
	d_filter = 5.e-3  # just by seeing a Giampaolo MM filter

	# absorptive loss of the mirror
	mirror_abs = 0.002  # from ref[3] of Charlies document
        mirror_rms = 0.2e-6 # mirror roughness achieved by planck
        epsilon_zero = 8.854e-12 #[F/m]
        mu_zero = 1.257e-6 #[H/m] 
        sigma_c = 36.9e6 #[S/m] Al conductivity
     
        

num_band = len(freq_GHz_arr)

if option_proposal == 'US_MO':
	# sky, HWP, aperture, lens1, lens2, filter, lenslet, det
	emiss_arr = [1., 0.03, 1., 0.1, 0.1, 0.1, 1., 1.]
	# sky, HWP-AR, aperture, baffle?, filter, lenslet, detector
	eff_arr = [1., 0.93, 1., 0.9, 0.9, 0.9, 0.98, 0.809]

if option_proposal == 'JAXA_MDR':
	# sky, HWP, aperture, mirror1, mirror2, filter, lenslet
	emiss_arr = [1., .03, 1.0, 0.005, 0.005, 0.1, 1., 1.000]
	# sky, HWP-AR, aperture, baffle?, filter, lenslet, detector
	eff_arr =   [1., .98, 1.0, 1.000, 1.000, 0.98, 0.809]


Tmis_year = 3.
Tmis_sec = Tmis_year*365.*24.*3600.
fsky = 1.

#uKarcmin_total_margin = 0.608
uKarcmin_total_margin = 1.

# detector specifications
#wafer_size = 117.e-3

T_bath = 0.10
T_bath_nominal = 0.1
current_noise = 7.e-12
X_fac = 2.5
#n_tc = 1.
n_tc = 3.	

#optimal_bolotemp_ratio = 2.14
optimal_bolotemp_ratio = 1.71 

	#n_tc = 3.
	#optimal_bolotemp_ratio = 1.705

	# Optics specifications
	#halfangle_edge_degs = 7.9  # F/#=3.6
#T_mir_arr = [4.]
#hwp_temperature_arr = np.array([5.0,7.5,10.,12.5,15.,17.5,20.,25.,30.])

T_mir1 = 5.
T_mir2 = 5.
T_ape = 2.
T_hwp = 5.
T_1K = 2.
T_lenslet = 0.1

T_mir1_nominal = 5.
T_mir2_nominal = 5.
T_ape_nominal = 2.
T_hwp_nominal = 5.
T_1K_nominal = 2.
T_lenslet_nominal = 0.1

# Toki's and Charlie's eff
def aperture_ext(freq_c, option='Name'):
	if int(freq_c*1e-9) == 60: i=0
	if int(freq_c*1e-9) == 78: i=1
	if int(freq_c*1e-9) == 100: i=2
	if int(freq_c*1e-9) == 140: i=3
	if int(freq_c*1e-9) == 195: i=4
	if int(freq_c*1e-9) == 280: i=5

	if option == 'Toki':
		Tedge_dB = np.array([3.17, 5.36, 8.81, 7.67, 14.9, 30.7])
		apt_eff = np.array([0.5181, 0.7088, 0.8684, 0.8291, 0.968, 0.999])

	return 0., Tedge_dB[i], apt_eff[i]
