import numpy as np
import pylab as py

pi = np.pi
c = 3.e8
radeg = (180./pi)
halfangle_edge_deg = 8.

def aperture(dpix,nu,halfangle_edge_deg):
    sigma0_rad = 1./2.*c/nu/pi*2.6/dpix
    Tedge_dB = - 10.* np.log10( np.exp(-(halfangle_edge_deg/radeg)**2/2. * (1./sigma0_rad)**2) )
    Pap = 1.-np.exp(-2.*0.115*Tedge_dB)
    return sigma0_rad, Tedge_dB, Pap

num = 100

#nu = 140.e9
dpix = np.arange(num)/float(num)*20.e-3 + 5.e-3

for nu in [60.,78.,100.,140.,195.,280.]:
    print nu
    sigma0_rad, Tedge_dB, Pap = aperture(dpix,nu*1.e9)
    
    py.subplot(311)
    py.plot(dpix*1e3, sigma0_rad*radeg)

    py.subplot(312)
    py.plot(dpix*1e3, Tedge_dB)
    py.ylim([0,30])

    py.subplot(313)
    py.plot(dpix*1e3, Pap)
    py.ylim([0,1.1])

py.subplot(311)
py.ylabel('Feed $\\sigma_0$ [degs]')

py.subplot(312)
py.ylabel('Edge taper [dB]')

py.subplot(313)
py.ylabel('Ap. eff.')
py.xlabel('Pixel diameter[mm]')
py.show()
