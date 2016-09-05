import numpy as np
import pylab as py
import lib_mappingspeed as lib_ms
import global_par as g
import os
import sys

pi = np.pi

freq_GHz = g.freq_GHz_arr
Pmax = 0.
dir_out = '/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20150701_LB_Sensitivity/data/tmp/'
bandwidth = g.bandwidth
margin = g.uKarcmin_total_margin
num_band = len(freq_GHz)

NET_out = []
NEP_out = []
w_1_out = []
Npix_out = []
Gave = []
Pload = []

for i in range(num_band):

	output = lib_ms.mapping_speed( freq_GHz[i], Pmax, dir_out, \
		bandwidth[i], g.wafer_size, g.aperture_diameter_mm, \
		g.emiss, g.eff, \
		g.T_bath, g.T_mir_arr[0], g.T_ape, g.T_hwp, g.T_1K, \
		g.T_bath_nominal, g.T_mir_nominal, g.T_ape_nominal, g.T_hwp_nominal, g.T_1K_nominal, g.halfangle_edge_degs)

	d_lens = output[10]
	if i<=2: ind = np.where((d_lens*1e3 > g.pix_low*(0.999)) & (d_lens*1e3 < g.pix_low*(1.001))) 
	if i>2: ind = np.where((d_lens*1e3 > g.pix_high*(0.999)) & (d_lens*1e3 < g.pix_high*(1.001))) 

	print d_lens[ind[0]]*1e3
	NET = output[32,ind[0]]
	NEP = output[26,ind[0]]
	Npix = np.int_(output[9,ind[0]]*g.wafer_num[i])

	if i>2: Npix = 370/2

	w_1 = 10800./np.pi * np.sqrt(8.*pi*NET**2/(g.Tmis_sec*Npix*2.))

	NEP_out.append(NEP)
	NET_out.append(NET)
	Npix_out.append(Npix)
	w_1_out.append(w_1)
	Gave.append(output[38,ind[0]])
	Pload.append(output[17,ind[0]])

print ''
for i in range(num_band):
	print '%1.2f %1.2f %1.0d  %1.3f aW/rtHz %1.3f %1.3f %1.2f    %1.3f pW  %1.3f pW/K' % \
		(g.T_bath, g.T_bath_nominal, Npix_out[i]*2, NEP_out[i][0]*1e18, NET_out[i][0], w_1_out[i][0], w_1_out[i][0]/margin, Pload[i][0]*1e12, Gave[i][0]*1e12)
tmp = np.sum((np.array(w_1_out[2:4])/margin)**(-2))
print 1./np.sqrt(tmp)
tmp = np.sum((np.array(w_1_out[:])/margin)**(-2))
print 1./np.sqrt(tmp)
print ''




