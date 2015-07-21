import numpy as np
import pynbody
import matplotlib.pyplot as plt
import anybins
def sigmaSFR(s,h,nbins=[10000,5000,2000,1000,500,200,100,50,20,10,5,2,1],halos=[1],tmin=0,tmax=False,dtave=1e9):
	'''
	calculates the standard deviation of SFRs averaged over different time intervals
	
	-----inputs-----
	s = simulation. used for s.properties object
	h = halo pynbody object (h = s.halos())
	nbins = list of number of bins to choose. When combined with tmin and tmax, this gives a list of dt over which to bin the SF . Should monotonically DECREASE
	halos = list of halos you want to do this for
	tmin/max = defines a range of times over which to do the calculation. tmin default is 0 and tmax default is the time associated with simulation snapshot s
	dtave = 	
	-----outputs-----
	sigmaSFR = an array of size len(halos) X ntimes with the standard deviations
	deltat = an array of size len(halos) X ntimes with the time intervals used
	'''
	if not tmax: tmax = s.properties['time'].in_units('yr')
	if not tmin: tmin = 0
	print tmax, tmin
	nhalos = len(halos)
	ntimes = len(nbins)
	sigmaSFR = np.zeros((nhalos,ntimes))
	print "bin numbers that will be used: ", nbins
	trange = float(tmax-tmin)
	nbins = np.array(nbins)
	dt = trange/nbins
	for ii in range(len(halos)):
		print "calculating data for halo number ", halos[ii]
		nh = halos[ii]
		SFMean = h[nh].stars['mass'].in_units('Msol').sum()/(s.properties['time'].in_units('yr').max() - h[nh].stars['tform'].in_units('yr').min())
		#tinit = s.properties['time'].in_units('yr')
		print "Mean SFR:", SFMean
		for jj in range(ntimes):
			print "dt = ", dt[jj], "number of steps:", nbins[jj]
			SFR,binedges = np.histogram(h[nh].stars['tform'][(h[nh].stars['tform']>0)].in_units('yr'),weights=h[nh].stars['massform'][(h[nh].stars['tform']>0)].in_units('Msol')/dt[jj],range=(tmin,tmax),bins=nbins[jj])
			bcenter = binedges[np.arange(len(binedges-1))]+0.5*dt[jj]
			SFAve = np.zeros(nbins[jj])
			for kk in range(nbins[jj]):
				o, = np.where((h[nh].stars['tform'][(h[nh].stars['tform']>0)].in_units('yr') > (bcenter[kk] - 0.5*dtave))&(h[nh].stars['tform'][(h[nh].stars['tform']>0)].in_units('yr') <= (bcenter[kk] + 0.5*dtave)))
				SFAve[kk] = h[nh].stars['massform'][o].in_units('Msol').sum()/dtave
				if bcenter[kk]-0.5*dtave < 0:
					dtactual = bcenter[kk] + 0.5*dtdave
				if bcenter[kk]+0.5*dtave > s.properties['time'].in_units('yr'):
					dtactual = s.properties['time'].in_units('yr') - (bcenter[kk] - 0.5*dtave)
			sigmaSFR[ii,jj] = np.std(np.log10(SFR[(SFR>0)]))
			
			
	return sigmaSFR, dt

def sigDeltaSFR(s,h,dTrange = [1e7,2e9], ndT = 10, logdT=True, halos=[1],ntimes=10000,massform=None):
	'''
	'''
	sigmadSFR = np.zeros((len(halos),ndT))
        sigmadSFRNorm = np.zeros((len(halos),ndT))
	sigmadSFRNorm2 = np.zeros((len(halos),ndT))
        sigmadSFRNorm3 = np.zeros((len(halos),ndT))
	if len(dTrange)!=2:
		print "dTrange must be a tuple of the form [low,high]"
		return
	if logdT==True: dT = dTrange[0] * 10**(np.arange(ndT)*float(np.log10(dTrange[1])-np.log10(dTrange[0]))/float(ndT-1))
	if logdT==False: dT = dTrange[0] + np.arange(ndT)*float(dTrange[1]-dTrange[0])/float(ndT-1)
	print "time ranges used: ", dT
	for ii in range(len(halos)):
		print "calculating data for halo number ", halos[ii]
		nh = halos[ii]
		for jj in range(ndT):
			print "binsize ", dT[jj]
			low=0.5*dT[jj]#h[nh].stars['tform'].min().in_units('yr')+0.5*dT[jj]
			high=s.properties['time'].in_units('yr')-dT[jj]
			if low > high:
				print "dT too high given star formation duration! dT < 2/3 * (total duration of SF for halo)"
				break
			dSFR = np.zeros(ntimes)
			times = np.random.uniform(low,high,size=ntimes)
			binsLow = times - 0.5*dT[jj]
			binsHigh = times + 0.5*dT[jj]
			bins1 = np.array([binsLow,binsHigh]).T
			if massform==None: weights = h[nh].stars['massform'].in_units('Msol')/dT[jj]
			else: 
				weights = np.zeros(len(h[nh].stars))
				weights[:] = massform/dT[jj]
			SFR1 = anybins.histogram(h[nh].stars['tform'].in_units('yr'),bins1,weights=weights)
			bins2 = bins1+dT[jj]
			SFR2 = anybins.histogram(h[nh].stars['tform'].in_units('yr'),bins2,weights=weights)
			deltaSFR = SFR2 - SFR1
			if massform==None: meanSFR = np.sum(h[nh].stars['massform'].in_units('Msol')/(h[nh].stars['tform'].in_units('yr').max()-h[nh].stars['tform'].in_units('yr').min()))
			else: meanSFR = np.sum(massform*len(h[nh].stars)/(h[nh].stars['tform'].in_units('yr').max()-h[nh].stars['tform'].in_units('yr').min()))
			#deltaSFRNorm = deltaSFR/(4.330647e-13*len(h[nh].stars)/(h[nh].stars['tform'].in_units('yr').max()-h[nh].stars['tform'].in_units('yr').min()))
			deltaSFRNorm = deltaSFR/meanSFR
			deltaSFRNorm2 = deltaSFR/deltaSFR.mean()
			deltaSFRNorm3 = deltaSFR/deltaSFR.sum()
			print deltaSFR.mean()
			print bins1
			sigmadSFR[ii,jj] = np.std(deltaSFR)
			sigmadSFRNorm[ii,jj] = np.std(deltaSFR)/meanSFR
			sigmadSFRNorm2[ii,jj] = np.std(deltaSFRNorm2)
			sigmadSFRNorm3[ii,jj] = np.std(deltaSFRNorm3)

	return sigmadSFR, sigmadSFRNorm, sigmadSFRNorm2,sigmadSFRNorm3,dT


def SBfrac(s,h,dtrange=[1e7,1e8],SBth = [2.,3.,4.],dtbins=10,massform=None):
	if massform:
		mf = np.ones(len(h.stars))*massform
	else:
		mf = h.stars['massform'].in_units('Msol')
	T = h.stars['tform'].in_units('yr').max()-h.stars['tform'].in_units('yr').min()
	Mstar = mf.sum()
	print Mstar
	sfrave = Mstar/T
	print sfrave
	logdt = (np.log10(dtrange[1])-np.log10(dtrange[0]))/np.float(dtbins)
	dtarr = np.arange(dtbins)*logdt + np.log10(dtrange[0])
	SBf = np.zeros((dtbins,len(SBth)))
	print "beginning loop"
	for i in range(dtbins):
		print i
		nbins = np.int(h.stars['tform'].in_units('yr').max()/(10**dtarr[i]))
		print nbins
		trange = [h.stars['tform'].in_units('yr').max() - nbins*10**dtarr[i],h.stars['tform'].in_units('yr').max()]
		print trange
		sfr,bins = np.histogram(h.stars['tform'].in_units('yr'),range=trange,bins=nbins,weights = mf/(10**dtarr[i]))
		print sfr.max()
		for j in range(len(SBth)):
			print j
			burst, = np.where(sfr>sfrave*SBth[j])
			print burst,np.sum(sfr[burst])
			SBf[i,j] = np.sum(sfr[burst]*10**dtarr[i])/Mstar

	return SBf,dtarr









	
	
