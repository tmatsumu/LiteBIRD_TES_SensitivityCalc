import numpy as np
import pylab as py

pi = np.pi
c = 3.e8

def absorption(thickness, losstan, index, freq):
	return np.exp(-2.*pi*losstan*thickness*index*freq/c)

losstan = 1e-4
index = 3.4 
freq = np.linspace(20,500,100)*1e9

py.plot(freq*1e-9, absorption(15.e-3, losstan, index, freq),label='15mm')
py.plot(freq*1e-9, absorption(20.e-3, losstan, index, freq),label='20mm')
py.plot(freq*1e-9, absorption(25.e-3, losstan, index, freq),label='25mm')
py.plot(freq*1e-9, absorption(30.e-3, losstan, index, freq),label='30mm')

y_lim = [0.9,1.02]
py.plot([40,40], y_lim, '--')
py.plot([50,50], y_lim, '--')
py.plot([68,68], y_lim, '--')
py.plot([78,78], y_lim, '--')
py.plot([89,89], y_lim, '--')
py.plot([100,100], y_lim, '--')
py.plot([140,140], y_lim, '--')
py.plot([166,166], y_lim, '--')
py.plot([195,195], y_lim, '--')
py.plot([235,235], y_lim, '--')

py.plot([280,280], y_lim, '--')
py.plot([337,337], y_lim, '--')
py.plot([402,402], y_lim, '--')

py.xlabel('Frequency [GHz]')
py.ylabel('Absorption')
py.ylim(y_lim)
py.xlim([20,500])
py.semilogx()
py.grid()
py.legend(loc='best')
py.show()