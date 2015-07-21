"""

bhanalysis
=====

by MJT
"""
import numpy as np
import matplotlib.pyplot as plt
from pynbody.analysis import profile, angmom, halo
from pynbody import filt, units, config, array
import warnings
import math

#from matplotlib import rc
#rc('text',usetex=True)


def BH_Disp(halos,nhalos=10,massive=False,title=False):
	for i in range(10):
		cen = halo.shrink_sphere_center(halos[i+1])
		halos[i+1]['pos'] -= cen
		bhind = np.where(halos[i+1].stars['tform'] < 0)
		print np.size(bhind)
		if not np.size(bhind): continue
		bhdist = halos[i+1].stars[bhind]['r'].in_units('kpc')
		Mhalo = halos[i+1]['mass'].in_units('Msol').sum()
		if (not massive):
			for ii in range(np.size(bhind)):
				plt.plot([Mhalo],[bhdist[ii]],'k.')
		else:
			if massive:
				mass = halos[i+1].stars[bhind]['mass']
				MassDist = bhdist[np.where(mass == mass.max())]
				print MassDist
				plt.plot([Mhalo],[MassDist],'k.')
				if title:
					plt.title(title)
				else:
					plt.title("Most Massive BHs")
	plt.ylabel("Distance from Halo Center (kpc)")
	plt.xlabel("Halo Mass (Msun)")

def BH_Disp_v_BHMass(halos,nhalos=10):
	for i in range(10):
		cen = halo.shrink_sphere_center(halos[i+1])
                halos[i+1]['pos'] -= cen
                bhind = np.where(halos[i+1].stars['tform'] < 0)
		if not np.size(bhind): continue
                bhdist = halos[i+1].stars[bhind]['r'].in_units('kpc')
		massbh = halos[i+1].stars[bhind]['mass'].in_units('Msol')
		plt.plot(massbh,bhdist,'k.')
	plt.ylabel("Displacement from Halo Center (kpc)")
	plt.xlabel("BH Mass (Msun)")
	

def BHHM_relation(halos,nhalos=10,center=True,massive=False):
	counter = 1
	while counter<=nhalos:
		print counter
		bhind = np.where(halos[counter].stars['tform'] < 0)
		if not np.size(bhind): continue
		Mhalo = halos[counter]['mass'].in_units('Msol').sum()
		Mbh = halos[counter].stars[bhind]['mass'].in_units('Msol')
		if massive:
			MassBH = Mbh[np.where(Mbh == Mbh.max())]
			plt.plot([Mhalo],MassBH,'k.')
		if center:
			cen = halo.shrink_sphere_center(halos[counter])
			halos[counter]['pos'] -= cen
			bhdist = halos[counter].stars[bhind]['r'].in_units('kpc')
			MassBH = Mbh[np.where(bhdist == bhdist.min())]
			plt.plot([Mhalo],MassBH,'b.')
		counter += 1
	plt.ylabel("BH Mass (Msun)")
	plt.xlabel("Halo Mass (Msun)")




	
