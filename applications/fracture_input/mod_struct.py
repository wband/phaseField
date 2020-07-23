import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit

def func(x, a, b):
    return 1/np.sqrt(2*np.pi*b)*np.exp(-np.square((a-x)/b)/2)

nx = 256
sx = 40

np.random.seed()
eta = np.random.normal(0.0,2.0,size=(nx,nx))


qx = np.arange(0,nx, dtype=np.float64)
qx = np.where(qx <= nx//2, qx/sx, (nx-qx)/sx)
qx *= 2 * np.pi
qy = np.arange(0,nx//2 +1, dtype=np.float64)
qy*= 2*np.pi/sx

q2 = (qx**2).reshape(-1,1) + (qy**2).reshape(1,-1)

filt = np.ones_like(q2)
# q points are (2*pi)
# filters lambda < 1
# q_s = 2*np.pi
# filters lambda < 2
qs = np.pi  #tag

filt *= (q2 <= qs ** 2)

h_qs = np.fft.irfftn(np.fft.rfftn(eta) * filt)
inp = np.zeros((nx+1,nx+1))
inp[0:nx,0:nx] = np.maximum(h_qs+1.0, 0.1)
inp[nx,:]=inp[0,:]
inp[:,nx]=inp[:,0]

nbins=20
hist, binedges = np.histogram(inp, bins=nbins)
bins = (binedges[0:nbins]+binedges[1:nbins+1])/2
histdb = np.double(hist)/np.double(hist.sum())
histdb = histdb/integrate.simps(histdb,bins)

popt, pcov = curve_fit(func, bins, histdb)

inp.tofile('data_sigma{0:f}_mu{1:f}.dat'.format(popt[1], popt[0]))
