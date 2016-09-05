import numpy as np
import pylab as py
import lib_m as lib_m
import lib_foreground as lib_f
import lib_Clmanip as libcl
import global_par as g
import func_foreground as func
import matplotlib.cm as cm
import sys

'''
combining 
	main_foreground_6band2component_survey4.py
	and
	mapping_speed_summary_plot.py

	optimization for 1, 2, 3, 4 color/pix
'''

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   define the pars
dir_data = '/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20150701_LB_Sensitivity/data/100mK_Dapt400mm/data/'
option_run = '1colorperpix'
#option_run = '2colorperpix'
#option_run = '3colorperpix'
filename_out = dir_data + '/monte_'+option_run
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   basic input parameters
pi = np.pi
margin = g.uKarcmin_total_margin
obs_time_sec = g.Tmis_sec

fsky = 0.53

nu_arr = g.freq_GHz_arr

A_tot = pi* 0.45**2/4.
N_tot_max_in = 1000
N_tot_min_in = 100
N_tot_max = 2500
N_tot_min = 1500

lmax = 20
num_iter = 100

option_plot = ''
option_unit = 'thermo'

num_band = g.num_band

r_in = g.r_in
nu1 = nu_arr[0]*1.e9
nu2 = nu_arr[1]*1.e9
nu3 = nu_arr[2]*1.e9
nu4 = nu_arr[3]*1.e9
nu5 = nu_arr[4]*1.e9
nu6 = nu_arr[5]*1.e9

filename1 = dir_data+'/mappingspeed_60GHz_F_T4.0K.npy'
filename2 = dir_data+'/mappingspeed_78GHz_F_T4.0K.npy'
filename3 = dir_data+'/mappingspeed_100GHz_F_T4.0K.npy'
filename4 = dir_data+'/mappingspeed_140GHz_F_T4.0K.npy'
filename5 = dir_data+'/mappingspeed_195GHz_F_T4.0K.npy'
filename6 = dir_data+'/mappingspeed_280GHz_F_T4.0K.npy'
out1 = np.load(filename1)
out2 = np.load(filename2)
out3 = np.load(filename3)
out4 = np.load(filename4)
out5 = np.load(filename5)
out6 = np.load(filename6)


num_monte = 10000

A_fp_arr = []
Npix_sum = []
uKarcmin_tot = []
r_wm_i_arr_m = []
r_wm_i_arr_s = []
uKarcmin1_arr = []
uKarcmin2_arr = []
uKarcmin3_arr = []
uKarcmin4_arr = []
uKarcmin5_arr = []
uKarcmin6_arr = []
Npix1_arr = []
Npix2_arr = []
Npix3_arr = []
Npix4_arr = []
Npix5_arr = []
Npix6_arr = []
dpix1_arr = []
dpix2_arr = []
dpix3_arr = []
dpix4_arr = []
dpix5_arr = []
dpix6_arr = []

for i in range(num_monte):
	if (i%1000) == 0: print i
#	print ''; print ''
	NET_to_w_arr = np.zeros(6)
	uKarcmin_arr = np.zeros(6)

	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if option_run == '1colorperpix':
		dpix_arr_ = (np.random.uniform(2.,25.,6))*1.e-3
		dpix_arr = np.array([dpix_arr_[0],dpix_arr_[1],dpix_arr_[2],dpix_arr_[3],dpix_arr_[4],dpix_arr_[5]])
		Npix_arr_ = np.random.uniform(N_tot_min_in,N_tot_max_in,6)
		Npix_arr = np.array([Npix_arr_[0],Npix_arr_[1],Npix_arr_[2],Npix_arr_[3],Npix_arr_[4],Npix_arr_[5]])
		for j in range(num_band):
			if j==0: out = out1
			if j==1: out = out2
			if j==2: out = out3
			if j==3: out = out4
			if j==4: out = out5
			if j==5: out = out6
			d_pix = out[10]
			del_dpix = d_pix[1]-d_pix[0]
			ind = np.where( (d_pix>dpix_arr[j]-del_dpix*0.5) & (d_pix<dpix_arr[j]+del_dpix*0.5) )
		#		Npix_arr[j] = np.float_(out[9,ind[0]])
			NET_to_w_arr[j] = out[32,ind[0]]
			uKarcmin_arr[j] = 10800./pi*np.sqrt(8.*pi*NET_to_w_arr[j]**2/(Npix_arr[j]*obs_time_sec))/margin

	if option_run == '2colorperpix':
		dpix_arr_ = (np.random.uniform(2.,25.,3))*1.e-3
		dpix_arr = np.array([dpix_arr_[0],dpix_arr_[0],dpix_arr_[1],dpix_arr_[1],dpix_arr_[2],dpix_arr_[2]])
		Npix_arr_ = np.random.uniform(N_tot_min_in,N_tot_max_in,3)
		Npix_arr = np.array([Npix_arr_[0],Npix_arr_[0],Npix_arr_[1],Npix_arr_[1],Npix_arr_[2],Npix_arr_[2]])
		for j in range(num_band):
			if j==0: out = out1
			if j==1: out = out2
			if j==2: out = out3
			if j==3: out = out4
			if j==4: out = out5
			if j==5: out = out6
			d_pix = out[10]
			del_dpix = d_pix[1]-d_pix[0]
			ind = np.where( (d_pix>dpix_arr[j]-del_dpix*0.5) & (d_pix<dpix_arr[j]+del_dpix*0.5) )
		#		Npix_arr[j] = np.float_(out[9,ind[0]])
			NET_to_w_arr[j] = out[32,ind[0]]
			uKarcmin_arr[j] = 10800./pi*np.sqrt(8.*pi*NET_to_w_arr[j]**2/(Npix_arr[j]*obs_time_sec))/margin

	if option_run == '3colorperpix':
		dpix_arr_ = (np.random.uniform(2.,25.,2))*1.e-3
		dpix_arr = np.array([dpix_arr_[0],dpix_arr_[0],dpix_arr_[0],dpix_arr_[1],dpix_arr_[1],dpix_arr_[1]])
		Npix_arr_ = np.random.uniform(N_tot_min_in,N_tot_max_in,2)
		Npix_arr = np.array([Npix_arr_[0],Npix_arr_[0],Npix_arr_[0],Npix_arr_[1],Npix_arr_[1],Npix_arr_[1]])
		for j in range(num_band):
			if j==0: out = out1
			if j==1: out = out2
			if j==2: out = out3
			if j==3: out = out4
			if j==4: out = out5
			if j==5: out = out6
			d_pix = out[10]
			del_dpix = d_pix[1]-d_pix[0]
			ind = np.where( (d_pix>dpix_arr[j]-del_dpix*0.5) & (d_pix<dpix_arr[j]+del_dpix*0.5) )
		#		Npix_arr[j] = np.float_(out[9,ind[0]])
			NET_to_w_arr[j] = out[32,ind[0]]
			uKarcmin_arr[j] = 10800./pi*np.sqrt(8.*pi*NET_to_w_arr[j]**2/(Npix_arr[j]*obs_time_sec))/margin
	out = 0
	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	Npix_sum_ = np.sum(Npix_arr)
	if ((Npix_sum_ > N_tot_max) | (Npix_sum_ < N_tot_min)): continue #print 'out count', Npix_sum_; continue

	A_fp_arr_ = np.sqrt(3.)/2. * np.sum(Npix_arr*dpix_arr**2)
	if (A_fp_arr_ > A_tot): continue #print 'out area'; continue

	uKarcmin1 = uKarcmin_arr[0]
	uKarcmin2 = uKarcmin_arr[1]
	uKarcmin3 = uKarcmin_arr[2]
	uKarcmin4 = uKarcmin_arr[3]
	uKarcmin5 = uKarcmin_arr[4]
	uKarcmin6 = uKarcmin_arr[5]

	FWHM1 = FWHM2 = FWHM3 = FWHM4 = FWHM5 = FWHM6 = 0.

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#   
	r_wm_i, alpha12_wm_i, alpha65_wm_i, alpha12in_wm_i, alpha65in_wm_i = func.func_foreground_6band2component_2dustcomp_v1( r_in,
																		nu1, nu2, nu3, nu4, nu5, nu6, 
																		uKarcmin1, uKarcmin2, uKarcmin3, uKarcmin4, uKarcmin5, uKarcmin6,
																		FWHM1, FWHM2, FWHM3, FWHM4, FWHM5, FWHM6,
																		lmax, num_iter,
																		option_unit,option_plot)

	uKarcmin_tot_ = 1./np.sqrt(1./(uKarcmin1**2) \
		+ 1./(uKarcmin2**2) \
		+ 1./(uKarcmin3**2) \
		+ 1./(uKarcmin4**2) \
		+ 1./(uKarcmin5**2) \
		+ 1./(uKarcmin6**2)) 

	if ((r_wm_i[1]/np.sqrt(fsky)) > 0.0005): continue #print 'out error', (r_wm_i[1]/np.sqrt(fsky)); continue

	A_fp_arr.append( A_fp_arr_ )
	Npix_sum.append(Npix_sum_)

	uKarcmin_tot.append(uKarcmin_tot_)

	uKarcmin1_arr.append(uKarcmin1)
	uKarcmin2_arr.append(uKarcmin2)
	uKarcmin3_arr.append(uKarcmin3)
	uKarcmin4_arr.append(uKarcmin4)
	uKarcmin5_arr.append(uKarcmin5)
	uKarcmin6_arr.append(uKarcmin6)

	Npix1_arr.append(Npix_arr[0])
	Npix2_arr.append(Npix_arr[1])
	Npix3_arr.append(Npix_arr[2])
	Npix4_arr.append(Npix_arr[3])
	Npix5_arr.append(Npix_arr[4])
	Npix6_arr.append(Npix_arr[5])

	dpix1_arr.append(dpix_arr[0])
	dpix2_arr.append(dpix_arr[1])
	dpix3_arr.append(dpix_arr[2])
	dpix4_arr.append(dpix_arr[3])
	dpix5_arr.append(dpix_arr[4])
	dpix6_arr.append(dpix_arr[5])

	r_wm_i_arr_m.append(r_wm_i[0] )
	r_wm_i_arr_s.append( r_wm_i[1]/np.sqrt(fsky) )
	print i+1, '/', num_monte

num = len(Npix_sum)
monte = range(num)

np.savez(filename_out+'.npz', \
	uKarcmin1_arr=uKarcmin1_arr, \
	uKarcmin2_arr=uKarcmin2_arr, \
	uKarcmin3_arr=uKarcmin3_arr, \
	uKarcmin4_arr=uKarcmin4_arr, \
	uKarcmin5_arr=uKarcmin5_arr, \
	uKarcmin6_arr=uKarcmin6_arr, \
	Npix1_arr=Npix1_arr, \
	Npix2_arr=Npix2_arr, \
	Npix3_arr=Npix3_arr, \
	Npix4_arr=Npix4_arr, \
	Npix5_arr=Npix5_arr, \
	Npix6_arr=Npix6_arr, \
	dpix1_arr=dpix1_arr, \
	dpix2_arr=dpix2_arr, \
	dpix3_arr=dpix3_arr, \
	dpix4_arr=dpix4_arr, \
	dpix5_arr=dpix5_arr, \
	dpix6_arr=dpix6_arr, \
	A_fp_arr=A_fp_arr, \
	Npix_sum=Npix_sum, \
	uKarcmin_tot=uKarcmin_tot, \
	r_wm_i_arr_m=r_wm_i_arr_m, \
	r_wm_i_arr_s=r_wm_i_arr_s, \
	fsky=fsky, \
	r_in=r_in, \
	A_tot=A_tot, \
	N_tot_max=N_tot_max, \
	nu_arr=nu_arr, 
	num_monte=num_monte, \
	monte=monte, \
	obs_time_sec=obs_time_sec, \
	help='uKarcmin1~6_arr, Npix1~6_arr, dpix1~6_arr, A_fp_arr, Npix_sum, uKarcmin_tot, r_wm_i_arr_m, r_wm_i_arr_s, fsky, r_in, A_tot, N_tot_max, nu_arr, num_monte, obs_time_sec, monte'
	)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.figure(figsize=(14,10))
py.subplot(421)
#py.plot( monte, Npix_sum, '.')
width = 0.35  
p1 = py.bar(monte, Npix1_arr, width, color='r')

Npix_accum = np.array(Npix1_arr)
p2 = py.bar(monte, Npix2_arr, width, color='g', bottom=Npix_accum)

Npix_accum = np.array(Npix2_arr) + Npix_accum
p3 = py.bar(monte, Npix3_arr, width, color='b', bottom=Npix_accum)

Npix_accum = np.array(Npix3_arr) + Npix_accum
p4 = py.bar(monte, Npix4_arr, width, color='c', bottom=Npix_accum)

Npix_accum = np.array(Npix4_arr) + Npix_accum
p5 = py.bar(monte, Npix5_arr, width, color='m', bottom=Npix_accum)

Npix_accum = np.array(Npix5_arr) + Npix_accum
p6 = py.bar(monte, Npix6_arr, width, color='y', bottom=Npix_accum)

py.ylabel('$N_{pix}$')
py.ylim([0,4000])
py.xlim([0,num+1])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(422)
py.plot( monte, A_fp_arr, 'o')
py.plot( monte, A_tot*np.ones(num), '-')
py.ylabel('$A_{fp}$ [m$^2$]')
py.ylim([0,A_tot*1.1])
py.xlim([0,num+1])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(423)
py.plot( monte, uKarcmin_tot, 'o')
py.ylabel('$w_{tot}^{-1}$ [$\mu$K.arcmin]')
py.ylim([0,10])
py.xlim([0,num+1])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(424)
py.errorbar( monte, r_wm_i_arr_m, r_wm_i_arr_s, fmt='o')
py.ylabel('$r$')
py.ylim([0,0.003])
py.xlim([0,num+1])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(425)
py.plot( monte, r_wm_i_arr_s, 'o')
py.ylabel('$\sigma_r$')
py.ylim([0,np.max(r_wm_i_arr_s)*1.1])
py.xlim([0,num+1])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(426)
#fig, ax = py.subplots(326)
#y = np.row_stack((uKarcmin1_arr, uKarcmin2_arr, uKarcmin3_arr, uKarcmin4_arr, uKarcmin5_arr, uKarcmin6_arr))
#py.stackplot( monte, uKarcmin1_arr, uKarcmin2_arr, uKarcmin3_arr, uKarcmin4_arr, uKarcmin5_arr, uKarcmin6_arr)
width = 0.35  
NETarr1 = np.array(uKarcmin1_arr/np.sqrt(Npix_sum))
p1 = py.bar(monte, NETarr1, width, color='r')

NETarr2 = np.array(uKarcmin2_arr/np.sqrt(Npix_sum))
NETarr_accum = NETarr1
p2 = py.bar(monte, NETarr2, width, color='g', bottom=NETarr_accum)

NETarr3 = np.array(uKarcmin3_arr/np.sqrt(Npix_sum))
NETarr_accum = NETarr_accum + NETarr2
p3 = py.bar(monte, NETarr3, width, color='b', bottom=NETarr_accum)

NETarr4 = np.array(uKarcmin4_arr/np.sqrt(Npix_sum))
NETarr_accum = NETarr_accum + NETarr3
p4 = py.bar(monte, NETarr4, width, color='c', bottom=NETarr_accum)

NETarr5 = np.array(uKarcmin5_arr/np.sqrt(Npix_sum))
NETarr_accum = NETarr_accum + NETarr4
p5 = py.bar(monte, NETarr5, width, color='m', bottom=NETarr_accum)

NETarr6 = np.array(uKarcmin6_arr/np.sqrt(Npix_sum))
NETarr_accum = NETarr_accum + NETarr5
p6 = py.bar(monte, NETarr6, width, color='y', bottom=NETarr_accum)

py.ylabel('$\sum_{\\nu}{NET_{arr}^{(\\nu)}}$')
py.xlim([0,num+1])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(427)
#fig, ax = py.subplots(326)
#y = np.row_stack((uKarcmin1_arr, uKarcmin2_arr, uKarcmin3_arr, uKarcmin4_arr, uKarcmin5_arr, uKarcmin6_arr))
#py.stackplot( monte, uKarcmin1_arr, uKarcmin2_arr, uKarcmin3_arr, uKarcmin4_arr, uKarcmin5_arr, uKarcmin6_arr)
width = 0.35  
dpix1_arr = np.array(dpix1_arr)*1e3
p1 = py.bar(monte, dpix1_arr, width, color='r')

dpix_arr_accum = np.array(dpix1_arr)
p2 = py.bar(monte, np.array(dpix2_arr)*1e3, width, color='g', bottom=dpix_arr_accum)

dpix_arr_accum = np.array(dpix2_arr)*1e3 + dpix_arr_accum
p3 = py.bar(monte, np.array(dpix3_arr)*1e3, width, color='b', bottom=dpix_arr_accum)

dpix_arr_accum = np.array(dpix3_arr)*1e3 + dpix_arr_accum
p4 = py.bar(monte, np.array(dpix4_arr)*1e3, width, color='c', bottom=dpix_arr_accum)

dpix_arr_accum = np.array(dpix4_arr)*1e3 + dpix_arr_accum
p5 = py.bar(monte, np.array(dpix5_arr)*1e3, width, color='m', bottom=dpix_arr_accum)

dpix_arr_accum = np.array(dpix5_arr)*1e3 + dpix_arr_accum
p6 = py.bar(monte, np.array(dpix6_arr)*1e3, width, color='y', bottom=dpix_arr_accum)

py.ylabel('$d_{pix}$ [mm]')
py.xlim([0,num+1])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.savefig(filename_out+'.png')

#py.show()

sys.exit()
