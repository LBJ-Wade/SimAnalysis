import numpy as np
import matplotlib.pyplot as plt
import pynbody
from pynbody.analysis import profile, angmom, halo, decomp
from pynbody import filt, units, config, array
from scipy import optimize as opt
import warnings
import math
import os
from pynbody.analysis import pkdgrav_cosmo as cosmo


def sigma(halo, sim, grps=[1],nbins=50, nangles=10,j_disk_min = 0.8, j_disk_max=1.1, E_cut = None, j_circ_from_r=False, vcen=None, log_interp=False, angmom_size="3 kpc"):
	sigma = np.zeros(np.size(grps))
	cen = np.zeros(np.size(grps))
	cnt = 0
	for grp in grps:
		print "calculating sigma for halo number ",grp 
		cen[cnt] = halo.shrink_sphere_center(h[grp])
		
