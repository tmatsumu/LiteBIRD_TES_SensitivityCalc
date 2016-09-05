import numpy as np
import pylab as py
import global_par as g
import sys
import os

pi = np.pi

dir_data = '/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20150701_LB_Sensitivity/data/100mK_Dapt400mm/data/'

num_mux = 6

#option_run = 'colorperpix_freeze1'
#option_run = 'colorperpix'
option_run = str(num_mux)+'colorperpix'

filename_out = dir_data + '/monte_'+option_run

out = np.load(filename_out+'.npz')

uKarcmin1_arr=out['uKarcmin1_arr']
uKarcmin2_arr=out['uKarcmin2_arr']
uKarcmin3_arr=out['uKarcmin3_arr']
uKarcmin4_arr=out['uKarcmin4_arr']
uKarcmin5_arr=out['uKarcmin5_arr']
uKarcmin6_arr=out['uKarcmin6_arr']
Npix1_arr=out['Npix1_arr']
Npix2_arr=out['Npix2_arr']
Npix3_arr=out['Npix3_arr']
Npix4_arr=out['Npix4_arr']
Npix5_arr=out['Npix5_arr']
Npix6_arr=out['Npix6_arr']
dpix1_arr=out['dpix1_arr']
dpix2_arr=out['dpix2_arr']
dpix3_arr=out['dpix3_arr']
dpix4_arr=out['dpix4_arr']
dpix5_arr=out['dpix5_arr']
dpix6_arr=out['dpix6_arr']
A_fp_arr=out['A_fp_arr']
Npix_sum=out['Npix_sum']
uKarcmin_tot=out['uKarcmin_tot']
r_wm_i_arr_m=out['r_wm_i_arr_m']
r_wm_i_arr_s=out['r_wm_i_arr_s']
fsky=out['fsky']
r_in=out['r_in']
A_tot=out['A_tot']
Ndet_tot_max=out['Ndet_tot_max']
nu_arr=out['nu_arr']
num_monte=out['num_monte']
monte=out['monte']
obs_time_sec=out['obs_time_sec']

if num_mux == 1: ind = np.where(r_wm_i_arr_s < 0.00054)
if num_mux == 2: ind = np.where(r_wm_i_arr_s < 0.000236)
#ind = np.where(r_wm_i_arr_s < 0.0001741)
if num_mux == 3: ind = np.where(r_wm_i_arr_s < 0.00016)
if num_mux == 6: ind = np.where(r_wm_i_arr_s < 0.00012)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

num_ = len(ind[0])
monte_ = np.linspace(1,num_,num_)

num = len(monte)

py.figure(1, figsize=(14,10))
py.subplot(421)
#py.plot( monte, Npix_sum, '.')
width = 0.35  
p1 = py.bar(monte_, Npix1_arr[ind[0]]*2, width, color='r')

Npix_accum = np.array(Npix1_arr)*2
p2 = py.bar(monte_, Npix2_arr[ind[0]]*2, width, color='g', bottom=Npix_accum[ind[0]])

Npix_accum = np.array(Npix2_arr)*2 + Npix_accum
p3 = py.bar(monte_, Npix3_arr[ind[0]]*2, width, color='b', bottom=Npix_accum[ind[0]])

Npix_accum = np.array(Npix3_arr)*2 + Npix_accum
p4 = py.bar(monte_, Npix4_arr[ind[0]]*2, width, color='c', bottom=Npix_accum[ind[0]])

Npix_accum = np.array(Npix4_arr)*2 + Npix_accum
p5 = py.bar(monte_, Npix5_arr[ind[0]]*2, width, color='m', bottom=Npix_accum[ind[0]])

Npix_accum = np.array(Npix5_arr)*2 + Npix_accum
p6 = py.bar(monte_, Npix6_arr[ind[0]]*2, width, color='y', bottom=Npix_accum[ind[0]])

Npix_accum = np.array(Npix6_arr)*2 + Npix_accum

py.ylabel('$N_{det}$')
py.ylim([0,Ndet_tot_max])
py.xlim([0,num_+1])
py.title(option_run)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fig = py.figure(1)
ax1 = fig.add_subplot(422)
ax2 = ax1.twinx()
#py.subplot(422)
ax1.plot( monte_, A_fp_arr[ind[0]], 'ob')
ax1.plot( monte_, A_tot*np.ones(num_), '-b')
ax1.set_ylabel('$A_{fp}$ [m$^2$]')
ax1.set_ylim([0,A_tot*1.1])
ax1.set_xlim([0,num_+1])

ax2.plot( monte_, np.sqrt(A_fp_arr[ind[0]]*4./pi), 'or')
ax2.plot( monte_, np.sqrt(A_tot*4./pi)*np.ones(num_), '-r')
ax2.set_ylabel('$D_{fp}$ [m]')
ax2.set_ylim([0,np.sqrt(A_tot*4./pi)*1.1])

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(423)
py.plot( monte_, uKarcmin_tot[ind[0]], 'o')
py.ylabel('$w_{tot}^{-1}$ [$\mu$K.arcmin]')
py.ylim([0,10])
py.xlim([0,num_+1])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(424)
py.errorbar( monte_, r_wm_i_arr_m[ind[0]], r_wm_i_arr_s[ind[0]], fmt='o')
py.ylabel('$r$')
py.ylim([-0.001,0.001])
py.xlim([0,num_+1])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(425)
py.plot( monte_, r_wm_i_arr_s[ind[0]], 'o')
py.ylabel('$\sigma_r$')
py.ylim([0,np.max(r_wm_i_arr_s[ind[0]])*1.1])
py.xlim([0,num_+1])
py.grid()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(426)
py.plot(monte_, uKarcmin1_arr[ind[0]], '.r')
py.plot(monte_, uKarcmin2_arr[ind[0]], '.g')
py.plot(monte_, uKarcmin3_arr[ind[0]], '.b')
py.plot(monte_, uKarcmin4_arr[ind[0]], '.c')
py.plot(monte_, uKarcmin5_arr[ind[0]], '.m')
py.plot(monte_, uKarcmin6_arr[ind[0]], '.y')
py.plot(monte_, uKarcmin_tot[ind[0]], 'ok')
#fig, ax = py.subplots(326)
#y = np.row_stack((uKarcmin1_arr, uKarcmin2_arr, uKarcmin3_arr, uKarcmin4_arr, uKarcmin5_arr, uKarcmin6_arr))
#py.stackplot( monte, uKarcmin1_arr, uKarcmin2_arr, uKarcmin3_arr, uKarcmin4_arr, uKarcmin5_arr, uKarcmin6_arr)
#width = 0.35  
#NETarr1 = np.array(uKarcmin1_arr/np.sqrt(Npix_sum))
#p1 = py.bar(monte, NETarr1, width, color='r')

#NETarr2 = np.array(uKarcmin2_arr/np.sqrt(Npix_sum))
#NETarr_accum = NETarr1
#p2 = py.bar(monte, NETarr2, width, color='g', bottom=NETarr_accum)

#NETarr3 = np.array(uKarcmin3_arr/np.sqrt(Npix_sum))
#NETarr_accum = NETarr_accum + NETarr2
#p3 = py.bar(monte, NETarr3, width, color='b', bottom=NETarr_accum)

#NETarr4 = np.array(uKarcmin4_arr/np.sqrt(Npix_sum))
#NETarr_accum = NETarr_accum + NETarr3
#p4 = py.bar(monte, NETarr4, width, color='c', bottom=NETarr_accum)

#NETarr5 = np.array(uKarcmin5_arr/np.sqrt(Npix_sum))
#NETarr_accum = NETarr_accum + NETarr4
#p5 = py.bar(monte, NETarr5, width, color='m', bottom=NETarr_accum)

#NETarr6 = np.array(uKarcmin6_arr/np.sqrt(Npix_sum))
#NETarr_accum = NETarr_accum + NETarr5
#p6 = py.bar(monte, NETarr6, width, color='y', bottom=NETarr_accum)

py.ylabel('$\sum_{\\nu}{NET_{arr}^{(\\nu)}}$')
py.xlim([0,num_+1])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(427)
if num_mux == 1:
	py.plot(monte_, dpix1_arr[ind[0]]*1e3, 'o')
	py.plot(monte_, dpix2_arr[ind[0]]*1e3, 'o')
	py.plot(monte_, dpix3_arr[ind[0]]*1e3, 'o')
	py.plot(monte_, dpix4_arr[ind[0]]*1e3, 'o')
	py.plot(monte_, dpix5_arr[ind[0]]*1e3, 'o')
	py.plot(monte_, dpix6_arr[ind[0]]*1e3, 'o')
if num_mux == 2:
	py.plot(monte_, dpix1_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix2_arr[ind[0]]*1e3, 'o')
	py.plot(monte_, dpix3_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix4_arr[ind[0]]*1e3, 'o')
	py.plot(monte_, dpix5_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix6_arr[ind[0]]*1e3, 'o')
if num_mux == 3:
	py.plot(monte_, dpix1_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix2_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix3_arr[ind[0]]*1e3, 'o')
	py.plot(monte_, dpix4_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix5_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix6_arr[ind[0]]*1e3, 'o')
if num_mux == 6:
	py.plot(monte_, dpix1_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix2_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix3_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix4_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix5_arr[ind[0]]*1e3, 'o')
#	py.plot(monte_, dpix6_arr[ind[0]]*1e3, 'o')
py.ylim([0,30])
py.ylabel('$d_{pix}$ [mm]')
py.xlim([0,num_+1])

#width = 0.35  
#dpix1_arr = np.array(dpix1_arr)*1e3
#p1 = py.bar(monte, dpix1_arr, width, color='r')
#
#dpix_arr_accum = np.array(dpix1_arr)
#p2 = py.bar(monte, np.array(dpix2_arr)*1e3, width, color='g', bottom=dpix_arr_accum)
#
#dpix_arr_accum = np.array(dpix2_arr)*1e3 + dpix_arr_accum
#p3 = py.bar(monte, np.array(dpix3_arr)*1e3, width, color='b', bottom=dpix_arr_accum)
#
#dpix_arr_accum = np.array(dpix3_arr)*1e3 + dpix_arr_accum
#p4 = py.bar(monte, np.array(dpix4_arr)*1e3, width, color='c', bottom=dpix_arr_accum)
#
#dpix_arr_accum = np.array(dpix4_arr)*1e3 + dpix_arr_accum
#p5 = py.bar(monte, np.array(dpix5_arr)*1e3, width, color='m', bottom=dpix_arr_accum)
#
#dpix_arr_accum = np.array(dpix5_arr)*1e3 + dpix_arr_accum
#p6 = py.bar(monte, np.array(dpix6_arr)*1e3, width, color='y', bottom=dpix_arr_accum)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fig = py.figure(1)
ax1 = fig.add_subplot(428)

ax1.plot(monte_, Npix_sum[ind[0]]*2,'ob')
ax1.plot(monte_, np.ones(num_)*Ndet_tot_max,'-b')
ax1.set_ylabel('$N_{pix}$')
ax1.set_ylim([0,Ndet_tot_max])
ax1.set_xlim([0,num_+1])

ax2 = ax1.twinx()
ax2.plot( monte_, r_wm_i_arr_s[ind[0]], 'hr')
ax2.plot( monte_, np.ones(num_)*0.00058, '-r')
ax2.set_ylabel('$\sigma_r$')
ax2.set_ylim([0,np.max(r_wm_i_arr_s[ind[0]])*1.1])

#ax1.plot(monte, Npix_sum*A_fp_arr*r_wm_i_arr_s,'ob')
#ax1.set_xlim([0,num+1])

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

idx = np.sort(Npix_accum[ind[0]])
#idx = np.sort(r_wm_i_arr_s[ind[0]])
ind_ = ind[0]
Npix_sort = []
idx_sort = []
j = 0
for i in sorted(enumerate(Npix_accum[ind[0]]), key=lambda x:x[1]):
#for i in sorted(enumerate(r_wm_i_arr_s[ind[0]]), key=lambda x:x[1]):
	Npix_sort.append(i[1])
	idx_sort.append(i[0])
	j += 1
	print j, i
	if j==20: break
print Npix_sort
ind = ind_[idx_sort]



#ind = np.where(r_wm_i_arr_s == np.min(r_wm_i_arr_s))
print ''
print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print '  Parameters for the minimum variance '
print ''
print 'monte#: ', monte[ind]
print ''
print 'Band: ', nu_arr
print ''
margin = g.uKarcmin_total_margin
print 'NET1: ', uKarcmin1_arr[ind]*(pi/10800.*margin) * np.sqrt((2.*Npix1_arr[ind]*obs_time_sec)/(8.*pi))
print 'NET2: ', uKarcmin2_arr[ind]*(pi/10800.*margin) * np.sqrt((2.*Npix2_arr[ind]*obs_time_sec)/(8.*pi))
print 'NET3: ', uKarcmin3_arr[ind]*(pi/10800.*margin) * np.sqrt((2.*Npix3_arr[ind]*obs_time_sec)/(8.*pi))
print 'NET4: ', uKarcmin4_arr[ind]*(pi/10800.*margin) * np.sqrt((2.*Npix4_arr[ind]*obs_time_sec)/(8.*pi))
print 'NET5: ', uKarcmin5_arr[ind]*(pi/10800.*margin) * np.sqrt((2.*Npix5_arr[ind]*obs_time_sec)/(8.*pi))
print 'NET6: ', uKarcmin6_arr[ind]*(pi/10800.*margin) * np.sqrt((2.*Npix6_arr[ind]*obs_time_sec)/(8.*pi))
print ''
print 'uKarcmin1: ', uKarcmin1_arr[ind]
print 'uKarcmin2: ', uKarcmin2_arr[ind]
print 'uKarcmin3: ', uKarcmin3_arr[ind]
print 'uKarcmin4: ', uKarcmin4_arr[ind]
print 'uKarcmin5: ', uKarcmin5_arr[ind]
print 'uKarcmin6: ', uKarcmin6_arr[ind]
print 'uKarcmin:  ', np.sqrt(1./(1./uKarcmin1_arr[ind]**2+1./uKarcmin2_arr[ind]**2+1./uKarcmin3_arr[ind]**2+1./uKarcmin4_arr[ind]**2+1./uKarcmin5_arr[ind]**2+1./uKarcmin6_arr[ind]**2))
print ''
print 'Npix1_arr: ', Npix1_arr[ind]
print 'Npix2_arr: ', Npix2_arr[ind]
print 'Npix3_arr: ', Npix3_arr[ind]
print 'Npix4_arr: ', Npix4_arr[ind]
print 'Npix5_arr: ', Npix5_arr[ind]
print 'Npix6_arr: ', Npix6_arr[ind]
print 'Npair total: ', Npix_sum[ind]
print 'Ndet total: ', Npix_sum[ind]*2
print ''
print 'dpix1_arr: ', dpix1_arr[ind]*1e3, 'mm'
print 'dpix2_arr: ', dpix2_arr[ind]*1e3, 'mm'
print 'dpix3_arr: ', dpix3_arr[ind]*1e3, 'mm'
print 'dpix4_arr: ', dpix4_arr[ind]*1e3, 'mm'
print 'dpix5_arr: ', dpix5_arr[ind]*1e3, 'mm'
print 'dpix6_arr: ', dpix6_arr[ind]*1e3, 'mm'
print ''
print 'delta r: ', r_wm_i_arr_s[ind]
print 'delta r_min: ', np.min(r_wm_i_arr_s)
print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print ''
#sys.exit()
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






num_ = len(ind)
monte_ = np.linspace(1,num_,num_)

py.figure(2, figsize=(7,8))
py.subplot(411)
if num_mux == 1:
	py.plot(monte_, dpix1_arr[ind]*1e3, 'or', markersize=9)
	py.plot(monte_, dpix2_arr[ind]*1e3, 'hg', markersize=9)
	py.plot(monte_, dpix3_arr[ind]*1e3, 'sb', markersize=9)
	py.plot(monte_, dpix4_arr[ind]*1e3, 'pc', markersize=9)
	py.plot(monte_, dpix5_arr[ind]*1e3, 'dm', markersize=9)
	py.plot(monte_, dpix6_arr[ind]*1e3, 'vy', markersize=9)
if num_mux == 2:
	py.plot(monte_, dpix1_arr[ind]*1e3, 'or', markersize=9)
#	py.plot(monte_, dpix2_arr[ind]*1e3, 'o', markersize=9)
	py.plot(monte_, dpix3_arr[ind]*1e3, 'hg', markersize=9)
#	py.plot(monte_, dpix4_arr[ind]*1e3, 'o', markersize=9)
	py.plot(monte_, dpix5_arr[ind]*1e3, 'sb', markersize=9)
#	py.plot(monte_, dpix6_arr[ind]*1e3, 'o', markersize=9)
if num_mux == 3:
	py.plot(monte_, dpix1_arr[ind]*1e3, 'or', markersize=9)
#	py.plot(monte_, dpix2_arr[ind]*1e3, 'o', markersize=9)
#	py.plot(monte_, dpix3_arr[ind]*1e3, 'o', markersize=9)
	py.plot(monte_, dpix4_arr[ind]*1e3, 'hg', markersize=9)
#	py.plot(monte_, dpix5_arr[ind]*1e3, 'o', markersize=9)
#	py.plot(monte_, dpix6_arr[ind]*1e3, 'o', markersize=9)
if num_mux == 6:
	py.plot(monte_, dpix1_arr[ind]*1e3, 'o')
	py.ylim([0,30])
py.xlim([0,num_+1])
py.grid()
py.ylabel('$d_{pix}$ [mm]', fontsize=13)
py.xticks(fontsize=13)
py.yticks(fontsize=13)
#py.title(option_run)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(412)
#py.plot( monte, Npix_sum, '.')
width = 0.4 
p1 = py.bar(monte_, Npix1_arr[ind]*2, width, color='r')

Npix_accum = np.array(Npix1_arr)*2
p2 = py.bar(monte_, Npix2_arr[ind]*2, width, color='g', bottom=Npix_accum[ind])

Npix_accum = np.array(Npix2_arr)*2 + Npix_accum
p3 = py.bar(monte_, Npix3_arr[ind]*2, width, color='b', bottom=Npix_accum[ind])

Npix_accum = np.array(Npix3_arr)*2 + Npix_accum
p4 = py.bar(monte_, Npix4_arr[ind]*2, width, color='c', bottom=Npix_accum[ind])

Npix_accum = np.array(Npix4_arr)*2 + Npix_accum
p5 = py.bar(monte_, Npix5_arr[ind]*2, width, color='m', bottom=Npix_accum[ind])

Npix_accum = np.array(Npix5_arr)*2 + Npix_accum
p6 = py.bar(monte_, Npix6_arr[ind]*2, width, color='y', bottom=Npix_accum[ind])

py.grid()
#py.ylim([0,N_tot_max])
py.ylim([0,4000])
py.xlim([0,num_+1])
py.ylabel('$N_{det}$', fontsize=13)
py.xticks(fontsize=13)
py.yticks(fontsize=13)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(413)
py.plot(monte_, uKarcmin1_arr[ind], 'or', markersize=9)
py.plot(monte_, uKarcmin2_arr[ind], 'hg', markersize=9)
py.plot(monte_, uKarcmin3_arr[ind], 'sb', markersize=9)
py.plot(monte_, uKarcmin4_arr[ind], 'pc', markersize=9)
py.plot(monte_, uKarcmin5_arr[ind], 'dm', markersize=9)
py.plot(monte_, uKarcmin6_arr[ind], 'vy', markersize=9)
py.plot(monte_, uKarcmin_tot[ind], 'ok', markersize=9)
py.ylim([0,25])
py.xlim([0,num_+1])
py.grid()
py.ylabel('$w_{tot}^{-1}$ [$\mu$K.arcmin]', fontsize=13)
py.xticks(fontsize=13)
py.yticks(fontsize=13)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#py.subplot(413)
#py.errorbar( monte_, r_wm_i_arr_m[ind[0]], r_wm_i_arr_s[ind[0]], fmt='o')
#py.ylabel('$r$')
#py.ylim([-0.001,0.001])
#py.xlim([0,num_+1])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
py.subplot(414)
py.plot( monte_, r_wm_i_arr_s[ind]*1e4, 'o', markersize=9)
py.ylabel('$\sigma_r$ [$\\times 10^{-4}$]')
#py.ylim([0,np.max(r_wm_i_arr_s[ind[0]]*1e4)*1.1])
py.ylim([0,8])
py.xlim([0,num_+1])
py.grid()
py.xticks(fontsize=13)
py.yticks(fontsize=13)
#py.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
py.xlabel('Candidates', fontsize=13)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#width = 0.35  
#dpix1_arr = np.array(dpix1_arr)*1e3
#p1 = py.bar(monte, dpix1_arr, width, color='r')
#
#dpix_arr_accum = np.array(dpix1_arr)
#p2 = py.bar(monte, np.array(dpix2_arr)*1e3, width, color='g', bottom=dpix_arr_accum)
#
#dpix_arr_accum = np.array(dpix2_arr)*1e3 + dpix_arr_accum
#p3 = py.bar(monte, np.array(dpix3_arr)*1e3, width, color='b', bottom=dpix_arr_accum)
#
#dpix_arr_accum = np.array(dpix3_arr)*1e3 + dpix_arr_accum
#p4 = py.bar(monte, np.array(dpix4_arr)*1e3, width, color='c', bottom=dpix_arr_accum)
#
#dpix_arr_accum = np.array(dpix4_arr)*1e3 + dpix_arr_accum
#p5 = py.bar(monte, np.array(dpix5_arr)*1e3, width, color='m', bottom=dpix_arr_accum)
#
#dpix_arr_accum = np.array(dpix5_arr)*1e3 + dpix_arr_accum
#p6 = py.bar(monte, np.array(dpix6_arr)*1e3, width, color='y', bottom=dpix_arr_accum)


py.show()