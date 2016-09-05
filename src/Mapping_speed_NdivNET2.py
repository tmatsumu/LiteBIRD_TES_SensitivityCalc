import numpy as np
import pylab as py
import exe_mappingspeed as exe_ms
#import os
import sys
import global_par as g

pi = np.pi
radeg = (180./pi)
c = 299792458.

print ""
print "Mapping_speed_NdivNET2"
print ""
dir_out = sys.argv[1]
filename_arr = exe_ms.mapping_speed( dir_out)

f1 = open("../data/"+g.option_proposal+"_par1.txt", "w")
f2 = open("../data/"+g.option_proposal+"_par2.txt", "w")
print ''
print 'Freq [GHz], D_lens [mm], apt., P_load [fW], NEP_ph, NEP_th, NEP_re, NEP_to, NET [uKrts], Npix, NETarr [uKrts], NETarr w/ m [uKrts]'

Ndet_final = []
num_band = len(filename_arr)
for i in range(num_band):
#	print filename_arr[i]
	output_tmp = np.load(filename_arr[i]+'.npy')
	output_arr = output_tmp.item()
	freq_GHz = output_arr['freq_GHz']
	keyname_arr = output_arr['keyname_arr']
	eff_arr = output_arr['eff_arr']
	emiss_arr = output_arr['emiss_arr']
	num_key = len(keyname_arr)
	output = output_arr['output']

	ind = np.where(np.array(keyname_arr) == 'D_lens')
	x = output[ind[0],:]
	del_x = x[0][1]-x[0][0]
	ind_key = np.where(np.array(keyname_arr) == 'P_load'); 	y_tmp = output[ind_key[0],:] 
	P_load = np.array(y_tmp[0])
	ind_key = np.where(np.array(keyname_arr) == 'D_lens'); 	y_tmp = output[ind_key[0],:] 
	D_lens = np.array(y_tmp[0])
	ind_key = np.where(np.array(keyname_arr) == 'apt_eff');  y_tmp = output[ind_key[0],:] 
	apt_eff = np.array(y_tmp[0])
	ind_key = np.where(np.array(keyname_arr) == 'NEP_ph_w');  y_tmp = output[ind_key[0],:] 
	NEP_ph_w = np.array(y_tmp[0])
	ind_key = np.where(np.array(keyname_arr) == 'NEP_th'); 	y_tmp = output[ind_key[0],:] 
	NEP_th = np.array(y_tmp[0])
	ind_key = np.where(np.array(keyname_arr) == 'NEP_readout');  y_tmp = output[ind_key[0],:] 
	NEP_readout = np.array(y_tmp[0])
	ind_key = np.where(np.array(keyname_arr) == 'NEP_to_w');  y_tmp = output[ind_key[0],:] 
	NEP_to_w = np.array(y_tmp[0])
	ind_key = np.where(np.array(keyname_arr) == 'NET_totalper_w'); 	y_tmp = output[ind_key[0],:] 
	NET_totalper_w = np.array(y_tmp[0])
	ind_key = np.where(np.array(keyname_arr) == 'Npix'); y_tmp = output[ind_key[0],:]
	Npix = np.array(y_tmp[0])
	ind_key = np.where(np.array(keyname_arr) == 'Vbias'); y_tmp = output[ind_key[0],:]
	Vbias = np.array(y_tmp[0])

#	if freq_GHz < g.freq_break:
	ind = np.where((x[0] > g.pixDmm[i]*1e-3-del_x/2.) & (x[0] < g.pixDmm[i]*1e-3+del_x/2.))
	ind = ind[0]
	print '  %1d %1.1f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1d %1.3f %1.3f %1.3f' % \
		(freq_GHz, \
		D_lens[ind[0]]*1e3, \
		apt_eff[ind[0]], \
		P_load[ind[0]]*1e12, \
		NEP_ph_w[ind[0]]*1e18, \
		NEP_th[ind[0]]*1e18, \
		NEP_readout[ind[0]]*1e18, \
		NEP_to_w[ind[0]]*1e18, \
		NET_totalper_w[ind[0]], \
		2.*Npix[ind[0]]*g.wafer_num[i]-g.stupid_offset[i], \
		NET_totalper_w[ind[0]]/np.sqrt(2.*Npix[ind[0]]*g.wafer_num[i]-g.stupid_offset[i]), \
		NET_totalper_w[ind[0]]/np.sqrt(2.*Npix[ind[0]]*g.wafer_num[i]-g.stupid_offset[i])*1.25, \
		NEP_readout[ind[0]]/Vbias[ind[0]]*1e12)

	f2.write( '%1d %1.1f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1d %1.3f %1.3f %1.3f \n' % \
		(freq_GHz, \
		D_lens[ind[0]]*1e3, \
		apt_eff[ind[0]], \
		P_load[ind[0]]*1e12, \
		NEP_ph_w[ind[0]]*1e18, \
		NEP_th[ind[0]]*1e18, \
		NEP_readout[ind[0]]*1e18, \
		NEP_to_w[ind[0]]*1e18, \
		NET_totalper_w[ind[0]], \
		2.*Npix[ind[0]]*g.wafer_num[i]-g.stupid_offset[i], \
		NET_totalper_w[ind[0]]/np.sqrt(2.*Npix[ind[0]]*g.wafer_num[i]-g.stupid_offset[i]), \
		NET_totalper_w[ind[0]]/np.sqrt(2.*Npix[ind[0]]*g.wafer_num[i]-g.stupid_offset[i])*1.25, \
		NEP_readout[ind[0]]/Vbias[ind[0]]*1e12))

	eff_arr[2] = apt_eff[ind[0]]
	num_tmp = len(emiss_arr)
	print ''
	print '  %s'  % ('freq_GHz, ref, eff_arr, emiss_arr')
	for j in range(num_tmp):
		print '  %1.1f %1.3f %1.3f %1.3f'  % (freq_GHz, g.ref_arr[j], eff_arr[j], emiss_arr[j])
		f1.write('  %1.1f %1.3f %1.3f %1.3f\n'  % (freq_GHz, g.ref_arr[j], eff_arr[j], emiss_arr[j]))
	f1.write('  \n')
	print ''

	Ndet_final.append(2.*Npix[ind[0]]*g.wafer_num[i])

	py.figure(1)
	ind = np.where(np.array(keyname_arr) == 'D_lens')
	x = output[ind[0],:]
	ind = np.where(np.array(keyname_arr) == 'P_load')
	y = output[ind[0],:]
	py.plot(x[0]*1e3,y[0]*1e12, linewidth=2, label=str(freq_GHz)+'GHz')
	py.xticks( color = 'k', size = 17)
	py.yticks( color = 'k', size = 17)
#	py.ylim([0,1.1])
	py.xlabel('Lens diameter [mm]', fontsize=17)
	py.ylabel('Total loading [fW]', fontsize=17)
	py.legend(loc='best')
	py.grid()

	py.figure(2)
	ind = np.where(np.array(keyname_arr) == 'D_lens')
	x = output[ind[0],:]
	ind = np.where(np.array(keyname_arr) == 'apt_eff')
	y = output[ind[0],:]
	py.plot(x[0]*1e3,y[0], linewidth=2, label=str(freq_GHz)+'GHz')
	py.xticks( color = 'k', size = 17)
	py.yticks( color = 'k', size = 17)
	py.ylim([0,1.1])
	py.xlabel('Lens diameter [mm]', fontsize=17)
	py.ylabel('Aperture efficiency', fontsize=17)
	py.legend(loc='best')
	py.grid()

	py.figure(3)
	ind = np.where(np.array(keyname_arr) == 'D_lens')
	x = output[ind[0],:]
	ind = np.where(np.array(keyname_arr) == 'NET_totalper_w')
	y = output[ind[0],:]
	py.plot(x[0]*1e3,y[0], linewidth=2, label=str(freq_GHz)+'GHz')
	py.xticks( color = 'k', size = 17)
	py.yticks( color = 'k', size = 17)
	py.ylim([0,200])
	py.xlabel('Lens diameter [mm]', fontsize=17)
	py.ylabel('$NET$ [$\mu$K$_{cmb}\sqrt{s}$]', fontsize=17)
	py.legend(loc=1)
	py.grid()

	py.figure(4)
	ind = np.where(np.array(keyname_arr) == 'D_lens')
	x = output[ind[0],:]
	ind = np.where(np.array(keyname_arr) == 'Npix')
	y = output[ind[0],:]
	py.plot(x[0]*1e3,2.*y[0]*g.wafer_num[i], linewidth=2, label=str(freq_GHz)+'GHz')
	py.xticks( color = 'k', size = 17)
	py.yticks( color = 'k', size = 17)
	py.ylim([0,150])
	py.xlabel('Lens diameter [mm]', fontsize=17)
	py.ylabel('$N_{det}$', fontsize=17)
	py.legend(loc='best')
	py.grid()

	py.figure(5)
	ind = np.where(np.array(keyname_arr) == 'D_lens')
	x = output[ind[0],:]
	ind = np.where(np.array(keyname_arr) == 'NET_totalper_w')
	y1 = output[ind[0],:]
	ind = np.where(np.array(keyname_arr) == 'Npix')
	y2 = output[ind[0],:]
	y = y1/np.sqrt(2.*y2*g.wafer_num[i])
	py.plot(x[0]*1e3,y[0], linewidth=2, label=str(freq_GHz)+'GHz')
	py.xticks( color = 'k', size = 17)
	py.yticks( color = 'k', size = 17)
	py.ylim([0,50])
	py.xlabel('Lens diameter [mm]', fontsize=17)
	py.ylabel('$NET_{arr}$ = $NET$/$\sqrt{N_{det}}$ [$\mu$K$_{cmb}\sqrt{s}$]', fontsize=17)
	py.legend(loc='best')

f1.close()
f2.close()

print Ndet_final
print np.sum(Ndet_final)

py.figure(1); py.grid()
py.figure(2); py.grid()
py.figure(3); py.grid()
py.figure(4); py.grid()
py.figure(5); py.grid()

#py.show()
sys.exit()


