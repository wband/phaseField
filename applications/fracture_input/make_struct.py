import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit

def func(x, a, b):
    return 1/np.sqrt(2*np.pi*b)*np.exp(-np.square((a-x)/b)/2)

nx = 256
sx = 40

eta = np.random.normal(0.0,1.0,size=(nx,nx))

inp = np.zeros((nx+1,nx+1))
inp[0:nx,0:nx] = eta
inp[nx,:]=inp[0,:]
inp[:,nx]=inp[:,0]

nbins=20
hist, binedges = np.histogram(inp, bins=nbins)
bins = (binedges[0:nbins]+binedges[1:nbins+1])/2
histdb = np.double(hist)/np.double(hist.sum())
histdb = histdb/integrate.simps(histdb,bins)

popt, pcov = curve_fit(func, bins, histdb)

inp.tofile('base_sigma{0:f}_mu{1:f}.dat'.format(popt[1], popt[0]))
