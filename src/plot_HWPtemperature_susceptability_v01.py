import numpy as np
import pylab as py
import sys
import fileinput

def read_txt13_senssum(filename):
	arr1=[]; arr2=[]; arr3=[]; arr4=[]; arr5=[]
	arr6=[]; arr7=[]; arr8=[]; arr9=[]; arr10=[]
	arr11=[]; arr12=[]; arr13=[]
	filelines = fileinput.input(filename)
	i=0
	for line in filelines:
		if i>=0:
			ar = line.split()
			arr1.append(int(ar[0]))
			arr2.append(float(ar[1]))
			arr3.append(float(ar[2]))
			arr4.append(float(ar[3]))
			arr5.append(float(ar[4]))
			arr6.append(float(ar[5]))
			arr7.append(float(ar[6]))
			arr8.append(float(ar[7]))
			arr9.append(float(ar[8]))
			arr10.append(int(ar[9]))
			arr11.append(float(ar[10]))
			arr12.append(float(ar[11]))
			arr13.append(float(ar[12]))
		i+=1
	return np.array(arr1),np.array(arr2),np.array(arr3),np.array(arr4),np.array(arr5), \
			np.array(arr6),np.array(arr7),np.array(arr8),np.array(arr9),np.array(arr10), \
			np.array(arr11),np.array(arr12),np.array(arr13)

dir_in = "/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20151223_LB_HFT_Sensitivity/"
dir_data = dir_in + "/data/Expected_HWPtemperature/data/HWPtemperature_susceptability/"
dir_out = dir_in + "/figures/Expected_HWPtemperature/HWPtemperature_susceptability/"

hwp_temperature_arr = np.array([5.0,7.5,10.,12.5,15.,17.5,20.,25.,30.])

num_i = len(hwp_temperature_arr)
num_band = 15

NETarr_sus40=[]
NETarr_sus50=[]
NETarr_sus60=[]
NETarr_sus68=[]
NETarr_sus78=[]
NETarr_sus89=[]
NETarr_sus100=[]
NETarr_sus119=[]
NETarr_sus140=[]
NETarr_sus166=[]
NETarr_sus195=[]
NETarr_sus235=[]
NETarr_sus280=[]
NETarr_sus337=[]
NETarr_sus402=[]

for i in range(num_i):
	tmp = str("%1.1f" % hwp_temperature_arr[i] )
#	Freq [GHz], D_lens [mm], apt., P_load [fW], NEP_ph, NEP_th, NEP_re, NEP_to, NET [uKrts], Npix, NETarr [uKrts], NETarr w/ m [uKrts]
	freq1, D_lens1, apt1, P_load1, NEP_ph1, NEP_th1, NEP_re1, NEP_to1, NET1, Npix1, NETarr1, NETarr_wm1, NEIarr1 \
		= read_txt13_senssum(dir_data+tmp+'K/US_MO_LFT_par2.txt')
	freq2, D_lens2, apt2, P_load2, NEP_ph2, NEP_th2, NEP_re2, NEP_to2, NET2, Npix2, NETarr2, NETarr_wm2, NEIarr2 \
		= read_txt13_senssum(dir_data+tmp+'K/US_MO_HFT_par2.txt')
	freq = np.hstack((freq1,freq2))
	D_lens = np.hstack((D_lens1,D_lens2))
	apt = np.hstack((apt1,apt2))
	P_load = np.hstack((P_load1,P_load2))
	NEP_ph = np.hstack((NEP_ph1,NEP_ph2))
	NEP_th = np.hstack((NEP_th1,NEP_th2))
	NEP_re = np.hstack((NEP_re1,NEP_re2))
	NEP_to = np.hstack((NEP_to1,NEP_to2))
	NET = np.hstack((NET1,NET2))
	Npix = np.hstack((Npix1,Npix2))
	NETarr = np.hstack((NETarr1,NETarr2))
	NETarr_wm = np.hstack((NETarr_wm1,NETarr_wm2))
	NEIarr = np.hstack((NEIarr1,NEIarr2))

	for j in range(num_band):
		if freq[j] == 40: NETarr_sus40.append(NETarr[j])
		if freq[j] == 50: NETarr_sus50.append(NETarr[j])
		if freq[j] == 60: NETarr_sus60.append(NETarr[j])
		if freq[j] == 68: NETarr_sus68.append(NETarr[j])
		if freq[j] == 78: NETarr_sus78.append(NETarr[j])
		if freq[j] == 89: NETarr_sus89.append(NETarr[j])
		if freq[j] == 100: NETarr_sus100.append(NETarr[j])
		if freq[j] == 119: NETarr_sus119.append(NETarr[j])
		if freq[j] == 140: NETarr_sus140.append(NETarr[j])
		if freq[j] == 166: NETarr_sus166.append(NETarr[j])
		if freq[j] == 195: NETarr_sus195.append(NETarr[j])
		if freq[j] == 235: NETarr_sus235.append(NETarr[j])
		if freq[j] == 280: NETarr_sus280.append(NETarr[j])
		if freq[j] == 337: NETarr_sus337.append(NETarr[j])
		if freq[j] == 402: NETarr_sus402.append(NETarr[j])

print len(hwp_temperature_arr), len(NETarr_sus40)
print hwp_temperature_arr, NETarr_sus40

py.figure()
py.subplot(211)
py.plot(hwp_temperature_arr,NETarr_sus40/NETarr_sus40[0],'-',label='40GHz')
py.plot(hwp_temperature_arr,NETarr_sus50/NETarr_sus50[0],'-',label='50GHz')
py.plot(hwp_temperature_arr,NETarr_sus60/NETarr_sus60[0],'-',label='60GHz')
py.plot(hwp_temperature_arr,NETarr_sus68/NETarr_sus68[0],'-',label='68GHz')
py.plot(hwp_temperature_arr,NETarr_sus78/NETarr_sus78[0],'-',label='78GHz')
py.plot(hwp_temperature_arr,NETarr_sus89/NETarr_sus89[0],'-',label='89GHz')
py.plot(hwp_temperature_arr,NETarr_sus100/NETarr_sus100[0],'-',label='100GHz')
py.plot(hwp_temperature_arr,NETarr_sus119/NETarr_sus119[0],'-',label='119GHz')
py.plot(hwp_temperature_arr,NETarr_sus140/NETarr_sus140[0],'-',label='140GHz')
py.plot(hwp_temperature_arr,NETarr_sus166/NETarr_sus166[0],'-',label='166GHz')
py.plot(hwp_temperature_arr,NETarr_sus195/NETarr_sus195[0],'-',label='195GHz')
py.plot(hwp_temperature_arr,NETarr_sus235/NETarr_sus235[0],'-',label='235GHz')
py.xlim([1,35])
py.ylim([0.9,2])
#py.xlabel('Temperature [K]',fontsize=17)
py.ylabel('$NET_{\\nu}(T)/NET(5K)$',fontsize=17)
py.xticks( color = 'k', size = 17)
py.yticks( color = 'k', size = 17)
#py.legend(loc='best')

py.subplot(212)
py.plot(hwp_temperature_arr,NETarr_sus280/NETarr_sus280[0],'-',label='280GHz')
py.plot(hwp_temperature_arr,NETarr_sus337/NETarr_sus337[0],'-',label='337GHz')
py.plot(hwp_temperature_arr,NETarr_sus402/NETarr_sus402[0],'-',label='402GHz')
py.xlim([1,35])
py.ylim([0.9,2])
py.xlabel('HWP Temperature [K]',fontsize=17)
py.ylabel('$NET_{\\nu}(T)/NET(5K)$',fontsize=17)
py.xticks( color = 'k', size = 17)
py.yticks( color = 'k', size = 17)
#py.legend(loc='best')

NET_combined = 1./np.sqrt(1./np.array(NETarr_sus40)**2  \
					    + 1./np.array(NETarr_sus50)**2 \
					    + 1./np.array(NETarr_sus60)**2 \
					    + 1./np.array(NETarr_sus68)**2 \
					    + 1./np.array(NETarr_sus78)**2 \
					    + 1./np.array(NETarr_sus89)**2 \
					    + 1./np.array(NETarr_sus100)**2 \
					    + 1./np.array(NETarr_sus119)**2 \
					    + 1./np.array(NETarr_sus140)**2 \
					    + 1./np.array(NETarr_sus166)**2 \
					    + 1./np.array(NETarr_sus195)**2 \
					    + 1./np.array(NETarr_sus235)**2 \
					    + 1./np.array(NETarr_sus280)**2 \
					    + 1./np.array(NETarr_sus337)**2 \
					    + 1./np.array(NETarr_sus402)**2 )

NET_cmb = 1./np.sqrt( 1./np.array(NETarr_sus78)**2 \
						+ 1./np.array(NETarr_sus89)**2 \
					    + 1./np.array(NETarr_sus100)**2 \
					    + 1./np.array(NETarr_sus119)**2 \
					    + 1./np.array(NETarr_sus140)**2 \
					    + 1./np.array(NETarr_sus166)**2 )

py.figure()
py.plot(hwp_temperature_arr,NET_combined/NET_combined[0],'-',label='Combined ch.')
py.plot(hwp_temperature_arr,NET_cmb/NET_cmb[0],'-',label='CMB ch.')
py.xlabel('Temperature [K]',fontsize=17)
py.ylabel('$NET(T)/NET(5K)$',fontsize=17)
py.xticks( color = 'k', size = 17)
py.yticks( color = 'k', size = 17)
py.xlim([1,35])
py.ylim([0.9,1.5])
py.legend(loc='best')

py.show()

sys.exit()














