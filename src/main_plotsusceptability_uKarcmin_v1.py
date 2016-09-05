import numpy as np
import pylab as py
import sys
import os
from scipy.optimize import curve_fit

'''
same as main_plotsusceptability_v1.py but the vertical axis is in uKarcmin_tot

'''

def cal_err(n,d,del_n,del_d):
	return np.sqrt((1./d)**2*del_n**2 + (n/d**2)**2*del_d**2)

def func_Tbath(x,par0,par1,par2):
	return par0+par1*np.exp(par2*(x-0.1))

def fit_Tbath(x,y,del_y,parin):
	par, cov = curve_fit(func_Tbath, x, y, sigma=del_y, p0=parin)
	print '->', par
	return par, cov

dir_in = '/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20150701_LB_Sensitivity/data/susceptability/rms20_LBbaseline/'
dir_out = '/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20150701_LB_Sensitivity/figures/'

#num_multichroic_arr = ['2', '3', '6']
#num_multichroic_arr = ['2', '3']
num_multichroic_arr = ['22', '33']
#num_multichroic_arr_name = ['2', '3']

#option_trial_arr = ['T_bath', 'T_mir', 'T_ape', 'T_hwp', 'T_1K', \
#				'T_bath_nominal', 'T_mir_nominal', 'T_ape_nominal', 'T_hwp_nominal', 'T_1K_nominal', \
#				'emiss_mir', 'eff_hwp', 'eff_filter', 'eff_lenslet', 'eff_det']
#option_trial_arr = ['T_1K','T_bath_nominal', 'T_mir_nominal', 'T_ape_nominal', 'T_hwp_nominal', 'T_1K_nominal', \
#				'emiss_mir', 'eff_hwp', 'eff_filter', 'eff_lenslet', 'eff_det']
#option_trial_arr = ['emiss_mir', 'eff_hwp', 'eff_filter', 'eff_lenslet', 'eff_det']
#option_trial_arr = ['eff_det']
option_trial_arr = ['T_bath', 'T_ape', 'T_hwp', 'T_bath_nominal', 'T_ape_nominal', 'T_hwp_nominal', 'eff_hwp', 'eff_filter', 'eff_lenslet', 'eff_det', 'margin']

fmt_c = ['r','g','b']

for option_trial in option_trial_arr:

	for num_multichroic in num_multichroic_arr:
		print num_multichroic, type(num_multichroic), type(dir_out), type(option_trial)
		filename_out = dir_in + '/multichroic'+num_multichroic+'_'+option_trial
		print filename_out+'.npz'
		output = np.load(filename_out+'.npz')

		r_arr=output['delta_r']
		syspar_arr=output['syspar_arr']
		uKarcmin_tot=output['uKarcmin_tot']
		nbdata = len(uKarcmin_tot)
		frac_err = np.ones(nbdata)
#		num_multichroic=output['num_multichroic']
#		option_trial=output['option_trial']

		if num_multichroic == '2': fmt_c_i = fmt_c[0]
		if num_multichroic == '3': fmt_c_i = fmt_c[1]
		if num_multichroic == '6': fmt_c_i = fmt_c[2]
		if num_multichroic == '22': fmt_c_i = fmt_c[0]; num_multichroic_name = '2'
		if num_multichroic == '33': fmt_c_i = fmt_c[1]; num_multichroic_name = '3'

		if option_trial == 'T_bath':
			py.subplot(211)
			print syspar_arr, len(syspar_arr)
			print uKarcmin_tot, len(uKarcmin_tot)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([0.08,0.16])
			py.ylim([2.8,3.8])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('$T_{bath}$', fontsize=18)

			x = np.linspace(0.08,0.22,100)
			py.subplot(212)
			parin = np.array([0.,1.,0.1])
#			parout, cov = fit_Tbath(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0],frac_err,parin)
#			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
#			py.plot(x,func_Tbath(x,parout[0],parout[1],parout[2]),'--'+fmt_c_i)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 1)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
#			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.plot(x,par[0]*x+par[1],'--'+fmt_c_i)
			py.xlim([0.08,0.16])
			py.ylim([0.98,1.1])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$T_{bath}$ [K] while the designed $T_{bath}$ is 0.1 [K]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

		if option_trial == 'T_ape':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([1,11])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('$T_{ape}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(1,12,100)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([1,11])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$T_{ape}$ [K] while the designed $T_{ape}$ is 4 [K]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

		if option_trial == 'T_hwp':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([1,11])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('$T_{hwp}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(1,12,100)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([1,11])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$T_{hwp}$ [K] while the designed $T_{hwp}$ is 4 [K]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

		if option_trial == 'T_1K':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([0.1,11])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('$T_{1K}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(0.,11,100)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*x+par[1],'--'+fmt_c_i)
#			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([0.1,11])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$T_{1K}$ [K] while the designed $T_{1K}$ is 1 [K]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

#++++++++++++++++++++++++++++++++++++++++++++++++

		if option_trial == 'T_bath_nominal':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([0.08,0.21])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('Baseline $T_{bath}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(0.04,0.21,100)
#			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 1)
#			par, cov = fit_Tbath(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0],frac_err,parin)
#			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)

			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)

#			py.plot(x,par[0]*x+par[1],'--'+fmt_c_i)
#			py.plot(x,func_Tbath(x,parout[0],parout[1],parout[2]),'--'+fmt_c_i)
			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([0.08,0.21])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$T_{bath}$ [K]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

		if option_trial == 'T_ape_nominal':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([1,11])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('Baseline $T_{ape}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(1,12,100)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([1,11])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$T_{ape}$ [K]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

		if option_trial == 'T_hwp_nominal':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([1,11])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('Baseline $T_{hwp}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(1,12,100)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([1,11])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$T_{hwp}$ [K]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

		if option_trial == 'T_1K_nominal':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([0.1,11])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('Baseline $T_{1K}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(0.,11,100)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*x+par[1],'--'+fmt_c_i)
#			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([0.1,11])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$T_{1K}$ [K]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

#++++++++++++++++++++++++++++++++++++++++++++++++

		if option_trial == 'eff_hwp':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([0.6,1.05])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('$\epsilon_{hwp}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(0.64,1.1,100)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([0.6,1.05])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$\epsilon_{hwp}$ while the designed $\epsilon_{hwp}$ is 0.98', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

		if option_trial == 'emiss_mir':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([0,0.055])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('$\epsilon_{mir}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(0,0.06,100)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*np.ones(len(x)),'--'+fmt_c_i)
#			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([0,0.055])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$\epsilon_{mir}$ while the designed $\epsilon_{mir}$ is 0.005', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

		if option_trial == 'eff_filter':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([0.55,1.05])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('$\epsilon_{filter}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(0.6,1.1,100)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([0.55,1.05])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$\epsilon_{filter}$ while the designed $\epsilon_{filter}$ is 0.95', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

		if option_trial == 'eff_det':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([0.35,0.95])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('$\epsilon_{det}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(0.3,1,100)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([0.35,0.95])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$\epsilon_{det}$ while the designed $\epsilon_{det}$ is 0.81', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

		if option_trial == 'eff_lenslet':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([0.55,1.05])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('$\epsilon_{lenslet}$', fontsize=18)

			py.subplot(212)
			x = np.linspace(0.6,1,100)
			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([0.55,1.05])
			py.ylabel('$w^{-1}$/$w^{-1}{(T_{bath}=0.1K)}$', fontsize=18)
			py.xlabel('$\epsilon_{lenslet}$ while the designed $\epsilon_{lenslet}$ is 0.99', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)

		if option_trial == 'margin':
			py.subplot(211)
			py.plot(syspar_arr, uKarcmin_tot[:], 'o'+fmt_c_i, label='$k=$'+num_multichroic_name)
			py.xlim([0.15,1.05])
			py.ylabel('$w^{-1}$ [$\mu$K.arcmin]', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.legend(loc='best')
			py.title('Scale factor', fontsize=18)

			py.subplot(212)
#			x = np.linspace(0.15,1.05,100)
#			par = np.polyfit(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 2)

			x = np.linspace(0.5,1.0,100)
			ind = np.where((syspar_arr<1.05) & (syspar_arr>0.5))
			par = np.polyfit(syspar_arr[ind[0]], uKarcmin_tot[ind[0]]/uKarcmin_tot[0], 2)
			print par[0], par[1], par[2]
			print '->', r_arr
			frac_err = cal_err(r_arr[0,:],r_arr[0,0],r_arr[1,:],r_arr[1,0])
			for i in range(len(syspar_arr)):
				print syspar_arr[i], uKarcmin_tot[i]/uKarcmin_tot[0], r_arr[0,i]*1e4, r_arr[0,i]/r_arr[0,0], frac_err[i]
			py.plot(syspar_arr, uKarcmin_tot[:]/uKarcmin_tot[0], 'o'+fmt_c_i)
			py.plot(x,par[0]*x**2+par[1]*x+par[2],'--'+fmt_c_i)
			py.xlim([0.15,1.05])
			py.ylabel('$w^{-1}$/$w^{-1}_0$', fontsize=18)
			py.xlabel('Scale factor', fontsize=18)
			py.xticks(fontsize=18)
			py.yticks(fontsize=18)
			py.grid()

	print option_trial
	py.savefig(dir_out+'/susceptability_uKarcmin_'+option_trial+'.eps')
	py.savefig(dir_out+'/susceptability_uKarcmin_'+option_trial+'.png')
	py.clf()

#	py.show()
sys.exit()

