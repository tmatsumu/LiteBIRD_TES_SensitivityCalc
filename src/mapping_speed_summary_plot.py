import numpy as np
import pylab as py
import global_par  as g
import lib_truncatedGaussian as lib_t
import sys

dir = '/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20150701_LB_Sensitivity/data/100mK_Dapt400mm/data/'

pi = np.pi

margin = g.uKarcmin_total_margin
obs_time_sec = g.Tmis_sec

nu_arr = g.freq_GHz_arr

Fnum = 3.5

py.figure(1,figsize=(8,14))
for nu_i in nu_arr:
	filename = dir+'mappingspeed_%1.0dGHz_F_T4.0K.npy' % (nu_i)
	out = np.load(filename)

	d_pix_mm = out[10]*1e3
	Npix = np.float_(out[9])
	NET_to_w = out[32]

	py.subplot(411)
	py.plot( d_pix_mm, NET_to_w)
	py.ylim([0,200])
	py.ylabel('$NET$ [$\mu$K$\sqrt{s}$]')

	if nu_i == nu_arr[0]:
		py.subplot(412)
		py.plot( d_pix_mm, Npix )
		py.ylim([0,100])
		py.ylabel('$N_{pix}$')

	py.subplot(413)
	py.plot( d_pix_mm, NET_to_w/np.sqrt(Npix) )
	py.ylim([0,50])
#	py.title('$NET_{arr} = NET/\sqrt{N_{pix}}$ [$\mu$K$\sqrt{s}$]')
	py.ylabel('$NET_{arr}$ [$\mu$K$\sqrt{s}$]')

	py.subplot(414)
	py.plot( d_pix_mm, 10800/pi*np.sqrt(8.*pi*NET_to_w**2/(Npix*obs_time_sec)) )
	py.ylim([0,50])
#	py.title('$\sqrt{8\pi} NET/\sqrt{N_{pix} t_{obs}}$ [$\mu$K.arcmin]')
	py.ylabel('$w^{-1}$ [$\mu$K.arcmin]')
	py.xlabel('Pixel diameter [mm]')


py.figure(2,figsize=(8,14))
for nu_i in nu_arr:
	filename = dir+'mappingspeed_%1.0dGHz_F_T4.0K.npy' % (nu_i)
	out = np.load(filename)

	d_pix_mm = out[10]*1e3
	Npix = np.float_(out[9])
	NET_to_w = out[32]

	py.subplot(411)
#	py.plot( d_pix_mm/Fnum, NET_to_w, label=str(nu_i)+'GHz')
	py.plot( d_pix_mm/Fnum, NET_to_w)
	py.ylim([0,200])
	py.ylabel('$NET$ [$\mu$K$\sqrt{s}$]', fontsize=16)
	py.xticks(fontsize=16)
	py.yticks(fontsize=16)
	py.legend(loc='best',prop={'size':12})

	if nu_i == nu_arr[0]:
		py.subplot(412)
		py.plot( d_pix_mm/Fnum, Npix )
		py.ylim([0,100])
		py.ylabel('$N_{pix}$', fontsize=16)
		py.xticks(fontsize=16)
		py.yticks(fontsize=16)

	py.subplot(413)
	py.plot( d_pix_mm/Fnum, NET_to_w/np.sqrt(Npix) )
	py.ylim([0,50])
#	py.title('$NET_{arr} = NET/\sqrt{N_{pix}}$ [$\mu$K$\sqrt{s}$]')
	py.ylabel('$NET_{arr}$ [$\mu$K$\sqrt{s}$]', fontsize=16)
	py.xticks(fontsize=16)
	py.yticks(fontsize=16)

	py.subplot(414)
	py.plot( d_pix_mm/Fnum, 10800/pi*np.sqrt(8.*pi*NET_to_w**2/(Npix*obs_time_sec)) )
	py.ylim([0,50])
#	py.title('$\sqrt{8\pi} NET/\sqrt{N_{pix} t_{obs}}$ [$\mu$K.arcmin]')
	py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=16)
	py.xlabel('Normalized pixel diameter, $d_{pix}/F\#$ [mm]', fontsize=16)
	py.xticks(fontsize=16)
	py.yticks(fontsize=16)


py.figure(3,figsize=(8,6))
for nu_i in nu_arr:
	filename = dir+'mappingspeed_%1.0dGHz_F_T4.0K.npy' % (nu_i)
	out = np.load(filename)

	d_pix_mm = out[10]*1e3
	Npix = np.float_(out[9])
	NET_to_w = out[32]
	apt_eff = out[5]

	Dapt_mm = 400.
	edge_dB = np.log(1.-apt_eff)/(-2.*0.115)

	num = len(edge_dB)
	num_ = 28
	j_arr = np.int_(np.linspace(10,num-1,20))
	theta_FHWM_arcmin = []
	theta_FHWM_edge_arcmin = []
	d_pix_ = []
	jj = 0
	for j in j_arr:
		filename_out = '/Users/tomotake_matsumura/Desktop/tmp/mappingspeed_%1.0dGHz_F_T4.0K%1d.png' % (nu_i, j)
		print j, num_
		if edge_dB[j] > 100: break
		d_pix_.append(d_pix_mm[j]/Fnum)
#		theta_FHWM_arcmin_, theta_FHWM_edge_arcmin_ = lib_t.lib_trancatedGaussianBeamShape(nu_i,Dapt_mm,-edge_dB[j],png_filenameout=filename_out)
		theta_FHWM_arcmin_, theta_FHWM_edge_arcmin_ = lib_t.lib_trancatedGaussianBeamShape(nu_i,Dapt_mm,-edge_dB[j],png_filenameout='')
		theta_FHWM_arcmin.append(theta_FHWM_arcmin_)
		theta_FHWM_edge_arcmin.append(theta_FHWM_edge_arcmin_)
		jj += 1

	py.subplot(211)
#	py.plot( d_pix_mm/Fnum, NET_to_w, label=str(nu_i)+'GHz')
	py.plot( d_pix_mm/Fnum, apt_eff, label=str(nu_i)+'GHz')
	py.ylim([0,1.1])
	py.ylabel('Aperture eff.', fontsize=16)
	py.xticks(fontsize=16)
	py.yticks(fontsize=16)
	py.legend(loc='best',prop={'size':12})

	py.subplot(212)
#	py.plot( d_pix_mm/Fnum, NET_to_w, label=str(nu_i)+'GHz')
#	if nu_i == 60: py.plot( d_pix_, theta_FHWM_edge_arcmin, 'or', label=str(nu_i)+'GHz', markersize=9)
#	if nu_i == 78: py.plot( d_pix_, theta_FHWM_edge_arcmin, 'hg', label=str(nu_i)+'GHz', markersize=9)
#	if nu_i == 100: py.plot( d_pix_, theta_FHWM_edge_arcmin, 'sb', label=str(nu_i)+'GHz', markersize=9)
#	if nu_i == 140: py.plot( d_pix_, theta_FHWM_edge_arcmin, 'pc', label=str(nu_i)+'GHz', markersize=9)
#	if nu_i == 195: py.plot( d_pix_, theta_FHWM_edge_arcmin, 'dm', label=str(nu_i)+'GHz', markersize=9)
#	if nu_i == 280: py.plot( d_pix_, theta_FHWM_edge_arcmin, 'vy', label=str(nu_i)+'GHz', markersize=9)
	if nu_i == 60: py.plot( d_pix_, theta_FHWM_edge_arcmin, label=str(nu_i)+'GHz', markersize=9)
	if nu_i == 78: py.plot( d_pix_, theta_FHWM_edge_arcmin, label=str(nu_i)+'GHz', markersize=9)
	if nu_i == 100: py.plot( d_pix_, theta_FHWM_edge_arcmin, label=str(nu_i)+'GHz', markersize=9)
	if nu_i == 140: py.plot( d_pix_, theta_FHWM_edge_arcmin, label=str(nu_i)+'GHz', markersize=9)
	if nu_i == 195: py.plot( d_pix_, theta_FHWM_edge_arcmin, label=str(nu_i)+'GHz', markersize=9)
	if nu_i == 280: py.plot( d_pix_, theta_FHWM_edge_arcmin,  label=str(nu_i)+'GHz', markersize=9)
#	py.ylim([0,200])
	py.ylabel('FWHM [arcmin]', fontsize=16)
	py.xticks(fontsize=16)
	py.yticks(fontsize=16)
	py.legend(loc='best',prop={'size':12})
	py.xlabel('Normalized pixel diameter, $d_{pix}/F\#$ [mm]', fontsize=16)


py.show()