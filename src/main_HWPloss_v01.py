import numpy as np
import pylab as py
import global_par as g

pi = np.pi
radeg = (180./pi)
c = 299792458.

freq_GHz = np.array(g.freq_GHz_arr)
freq_c = freq_GHz*1.e9
loss_HWP_LFT = 1. - g.ref_arr[1] - (1. - np.exp(-2.*pi*freq_c/c*g.d_lens*g.n_hwp*g.losstan_hwp))

py.plot(freq_GHz,loss_HWP_LFT)
py.ylim([0,1])

print np.mean(loss_HWP_LFT)
py.show()