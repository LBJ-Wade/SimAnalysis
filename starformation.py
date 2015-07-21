"""

starformation
=====

update to some of the functions in the stars.py file
by MJT
"""

import numpy as np
import matplotlib.pyplot as plt
import pynbody
from pynbody.analysis import profile, angmom, halo
from pynbody import filt, units, config, array
import warnings
import math
import anybins
import bhanalysis
plt.ion()
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
plt.rc('font', weight='medium')
plt.rc('axes', linewidth=2)
plt.rc('xtick.major',width=2)
plt.rc('ytick.major',width=2)
def partial_derivative(func, var=0, point=[]):
    args = point[:]
    def wraps(x):
        args[var] = x
        return func(*args)
    return scipy.misc.derivative(wraps, point[var], dx = 1e-8)

def sfh(sim,filename=None,massform=True,initstarmass=False,makeplot=True,
        subplot=False, trange=False, bins=100, binsize=False, zmin=False,overplot=False,**kwargs):

    '''
    star formation history

    **Optional keyword arguments:**

       *trange*: list, array, or tuple
         size(t_range) must be 2. Specifies the time range.

       *bins*: int
         number of bins to use for the SFH

       *label*: string
         label line for legend

	*zmin*: float
	 set min z to plot on second axes

	*overplot*: bool
	 set to True if you are plotting a line on an already existing plot


       *massform*: bool
         decides whether to use original star mass (massform) or final star mass

       *subplot*: subplot object
         where to plot SFH

       *legend*: boolean
         whether to draw a legend or not

    By default, sfh will use the formation mass of the star.  In tipsy, this will be
    taken from the starlog file.  Set massform=False if you want the final (observed)
    star formation history

    **Usage:**

    >>> import pynbody.plot as pp
    >>> pp.sfh(s,linestyle='dashed',color='k')


    '''

    if subplot:
        plt = subplot
    else:
	import matplotlib.pyplot as plt

    if 'nbins' in kwargs:
        bins=kwargs['nbins']
        del kwargs['nbins']

    if trange:
        assert len(trange) == 2
    else:
        trange = [0,sim.star['tform'].in_units("Gyr").max()]
    if binsize:
	bins = int((trange[1] - trange[0])/binsize)
	binnorm = 1./(1e9*binsize)
    else:
    	binnorm = 1e-9*(bins / (trange[1] - trange[0]))
    trangefilt = filt.And(filt.HighPass('tform',str(trange[0])+' Gyr'),
                          filt.LowPass('tform',str(trange[1])+' Gyr'))
    
    tforms = sim.star[trangefilt]['tform'].in_units('Gyr')
    if massform and not initstarmass:
        try:
	    weight = sim.star[trangefilt]['massform'].in_units('Msol') * binnorm
        except (KeyError, units.UnitsException) :
            warnings.warn("Could not load massform array -- falling back to current stellar masses", RuntimeWarning)
	    weight = sim.star[trangefilt]['mass'].in_units('Msol') * binnorm
    if initstarmass:
	    weight = np.zeros(np.size(tforms))
	    weight[:] = initstarmass*binnorm
    if not initstarmass and not massform:
        weight = sim.star[trangefilt]['mass'].in_units('Msol') * binnorm


    if not makeplot:
    	sfhist, thebins = np.histogram(tforms, weights=weight, range=trange,bins=bins)
    if makeplot:
	   sfhist, thebins, patches = plt.hist(tforms, weights=weight, range=trange,bins=bins,histtype='step',**kwargs)
 	   if not overplot:
 	   	if not subplot:
 		       	plt.ylim(0.0,1.2*np.max(sfhist))
			plt.xlim(trange)
        		plt.xlabel('Time [Gyr]',fontsize='large')
        		plt.ylabel('SFR [M$_\odot$ yr$^{-1}$]',fontsize='large')
    		else:
        		plt.set_ylim(0.0,1.2*np.max(sfhist))

	    # Make both axes have the same start and end point.
	   if subplot: x0,x1 = plt.get_xlim()
	   else: x0,x1 = plt.gca().get_xlim()
	   from pynbody.analysis import pkdgrav_cosmo as cosmo
	   c = cosmo.Cosmology(sim=sim)
	    
	   if not overplot:
	    	pz = plt.twiny()
	    	if not zmin:
			labelzs = np.arange(5,int(sim.properties['z'])-1,-1)
	    	else:
			labelzs = np.arange(10,int(sim.properties['z'])-1,-1)
		times = [13.7*c.Exp2Time(1.0 / (1+z))/c.Exp2Time(1) for z in labelzs]
	    	pz.set_xticks(times)
	    	pz.set_xticklabels([str(x) for x in labelzs])
	    	pz.set_xlim(x0, x1)
	    	pz.set_xlabel('z')
	
	   if (filename):
	        if config['verbose']: print "Saving "+filename
	        plt.savefig(filename)


    return array.SimArray(sfhist, "Msol yr**-1"), array.SimArray(thebins, "Gyr")



def genCSFRfit(z,z0,A,B,C):
	return  C/(10**(A*(z-z0)) + 10**(B*(z-z0)))

def CSFRFit(z,type='beh'):
	if type=='beh':
		#Behroozi 13
		z0 = 1.243
		C = 0.18
		A = -0.997
		B = 0.241
	if type=='hop':
		#Hopkins 06
		z0 = 0.840
                C = 0.143
                A = -1.311
                B = 0.085
	sigma = np.zeros(len(z))
	sigma[(z<=0.9)] = 0.13
	sigma[((0.9<z)&(z<=1.5))] = 0.17
	sigma[((1.5<z)&(z<=3))] = 0.19
	sigma[(3<z)] = 0.27
        return genCSFRfit(z,z0,A,B,C), sigma



def plotCSFRdata(Volume):
	#Duncan 14

	zz = np.array([4, 5, 6, 7])+1
	ss = np.array([-1.14,-1.33,-1.58,-1.78])
	plt.scatter(zz, 10**ss, color='k', marker='D', label='Duncan, K et al 2014', s=40)
	#Kistler 13:
	sfr = np.array([0.0653, 0.03, 0.041, 0.0276, 0.025])
	zhigh = [4.5, 5.5, 6.75, 8, 9.4]
	zhighplus = np.array([0.5, 0.5, 0.75, 0.5, 1.025])+1
	zhighminus = np.array([0.5, 0.5, 0.75, 0.5, 1.025])+1
	sfrplus = np.array([0.0653, 0.03, 0.0405, 0.0647, 0.058])
	sfrminus = np.array([0.0326, 0.015, 0.023, 0.023, 0.021])
	logsfrplus = np.log10(sfrplus + sfr) - np.log10(sfr)
	logsfrminus = np.abs(np.log10(sfr - sfrminus) - np.log10(sfr))
	plt.errorbar(zhigh, sfr, fmt='o',yerr=[sfrminus, sfrplus], xerr=[zhighminus, zhighplus], ls='None', linewidth=1.5, color='k', label='Kistler+ 13')
	#plt.scatter(zhigh, np.log10(sfr),marker='o', s=80, linewidth=0, color='orange')
	#Behroozi 13 and Hopkins 06 fitted relations
	zfits = np.arange(1,11,0.01)
	fitB,sigB = CSFRFit(zfits,type='beh')
	fitH,sigH = CSFRFit(zfits,type='hop')
	plt.plot(zfits+1,fitB,'k-',label='Behroozi+ 13')
	plt.plot(zfits+1,fitH,'k--',label='Hopkins+ 06')
	plt.fill_between(zfits+1,10**(np.log10(fitB)-sigB),10**(np.log10(fitB)+sigB),linewidth=1.5,facecolor='grey',alpha=0.2)
	plt.fill_between(zfits+1,10**(np.log10(fitH)-sigH),10**(np.log10(fitH)+sigH),linewidth=1.5,facecolor='grey',alpha=0.2)
	return
def cosmicSF(sim, sl, bins=100,Volume=50.**3,zrange=False,logbins=True,massform=True,initmass=False):
	from pynbody.analysis import cosmology
	if zrange:
     		assert len(zrange) == 2
    	else:
        	zrange = [0,20]
	if logbins==True:
	        logbinsize = (np.log10(zrange[1]+1)-np.log10(zrange[0]+1))/bins
	        logzplusonebins = np.log10(zrange[0]+1)+np.arange(bins+1)*logbinsize
	        zplusonebins = 10**logzplusonebins
	        zbins = zplusonebins-1
    	if not logbins:
        	zbinsize = np.float(zrange[1]-zrange[0])/bins
        	zbins = zrange[0] +np.arange(bins+1)*zbinsize
   	from pynbody.analysis import pkdgrav_cosmo as cosmo
    	c = cosmo.Cosmology(sim=sim)
	print zbins
    	timebins = np.zeros(bins+1)
    	for ii in range(bins+1):
       		timebins[ii] = 13.7-13.7*c.Exp2Time(1.0 / (1+zbins[ii]))/c.Exp2Time(1)
	if massform and not initmass:
		weight = sl['massform'].in_units('Msol')
	if initmass:
		weight = np.zeros(len(sl['tform']))
		weight[:] = initmass
	print timebins
	SF,binedges = np.histogram(13.7-sl['tform'].in_units('Gyr'),bins=timebins,weights=weight)
	binnorm = np.zeros(bins)
	for i in range(bins):
                binnorm[i] = 1e-9 / (timebins[i+1] - timebins[i])
	print SF, binnorm
	sfrdens = SF*binnorm/Volume
	return sfrdens, timebins, zbins
				
def plotCSFR(s,sl,bins=100,Volume=50.**3,massform=True,initmass=False,log=True,zrange=False,overplot=False,style='k-',label=None,plotdata=True,linewidth=4):
	sfrdens, timebins, zbins = cosmicSF(s,sl,bins=bins,Volume=Volume,massform=massform,initmass=initmass,zrange=zrange,logbins=log)
	plt.plot((zbins[0:-1]+0.5*(zbins[1:]-zbins[0:-1]))+1,sfrdens,style,label=label,linewidth=linewidth)
	if log: 
		plt.yscale('log',base=10)
		plt.xscale('log',base=10)
	if not overplot:
		plt.xticks([1,2,3,4,5,6,7,8,9,10],['0','1','2','3','4','5','6','7','8','9'])
	if not overplot and plotdata==True:
		plotCSFRdata(Volume)
	return
	
def plotSFHandBH(h,s,BHorbit,axarr=None,BHids=[],SFbins=100,trange=None,BHbins=100,tunits='Gyr',label=None,overplot=False,SFcolor='b',BHcolor=['k'],SFline='solid',BHline=['solid'],plotBHDetail=True,maxBH=False):
	if trange==None:
                trange = [0,s.properties['time'].in_units('Gyr')]
		trangeA = pynbody.array.SimArray(trange,tunits)
        dtyr = (trangeA.in_units('yr')[1] - trangeA.in_units('yr')[0])/SFbins
	if axarr != None and overplot==False: overplot=True
	if axarr == None: f,axarr = plt.subplots(2,sharex=True)
	if len(axarr) != 2: print "WARNING: more subplots than two were defined!"
	#SFH (top plot)
	print "getting SFH..."
	axarr[0].hist(h.stars['tform'][(h.stars['tform']>0)].in_units(tunits),bins=SFbins,range=trange,weights=h.stars['massform'][(h.stars['tform']>0)].in_units('Msol')/dtyr,color=SFcolor,histtype='step',linestyle=SFline,label=label,lw=2)
	print "getting BH accretion history for all BHs"

	cnt = 0	
	if maxBH==True and plotBHDetail==True:
		print "WARNING not obvious how to plot with both plotBHDetail and BHmax... just doing BHmax"
		plotBHDetail=False
	if maxBH == True: smaccList=[]
	for id in BHids:
		o, = np.where(BHorbit['iord']==id)
		print "BH ID ", id
		if len(o)==0:
			print "BH id", id, " not found in orbit file..."
			continue
		bho = BHorbit['data'][o[0]]
		smacc,times = bhanalysis.smoothAcc(bho,trange=trange,bins=BHbins,tunits=tunits)
		fa,= np.where(smacc>0)
		st = fa[0]
		fin = fa[-1]
		tfa = times[st]
		if cnt < len(BHids)-1: 
			if maxBH==False: axarr[1].plot(times[st:fin],smacc[st:fin]*0.1*3e10*3e10,color=BHcolor[cnt],linestyle=BHline[cnt],lw=2)
			else: smaccList.append(smacc)
		else: 
			if maxBH==False: axarr[1].plot(times[st:],smacc[st:]*0.1*3e10*3e10,color=BHcolor[cnt],linestyle=BHline[cnt],lw=2)
			else:
				smaccList.append(smacc)	
				smaccplot = np.array(smaccList).max(axis=0)
				fa,= np.where(smacc>0)
                		st = fa[0]
				axarr[1].plot(times[st:],smaccplot[st:]*0.1*3e10*3e10,color=BHcolor[cnt],linestyle=BHline[cnt],lw=2)
		if plotBHDetail==True:
			if cnt < len(BHids)-1: axarr[1].plot(bho['Time'][((bho['Time'].in_units(tunits)<times[-2])&(bho['Time'].in_units(tunits)>tfa))].in_units(tunits),bho['mdot'][((bho['Time'].in_units(tunits)<times[-2])&(bho['Time'].in_units(tunits)>tfa))].in_units('g s**-1')*0.1*3e10*3e10,color=BHcolor[cnt],alpha=0.25,linestyle='-')
			else: axarr[1].plot(bho['Time'][(bho['Time'].in_units(tunits)>tfa)].in_units(tunits),bho['mdot'][(bho['Time'].in_units(tunits)>tfa)].in_units('g s**-1')*0.1*3e10*3e10,color=BHcolor[cnt],alpha=0.25,linestyle='-')
		cnt += 1
	if not overplot:
		axarr[0].set_ylabel(r'SFR [M$_{\odot}$ yr$^{-1}$]',fontsize=30)
		plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
		axarr[1].set_yscale('log',base=10)
		axarr[1].set_ylabel(r'L$_{bol}$ [ergs/s]',fontsize=30)
		axarr[1].set_xlabel(r'Time (Gyr)',fontsize=30)
		plt.subplots_adjust(hspace=0)
	if overplot==False: return axarr,f
	else: return axarr
	
		


def sfdens(sim,Volume=50.**3, filename=None,massform=True,clear=True,legend=False,
        subplot=False, bins=100, label=False, zrange=False,overplot=False,logbins=True,histogram=False,pltloglog=True,**kwargs):

    '''
    star formation history

    **Optional keyword arguments:**

       *trange*: list, array, or tuple
         size(t_range) must be 2. Specifies the time range.

       *nbins*: int
         number of bins to use for the SFH

       *massform*: bool
         decides whether to use original star mass (massform) or final star mass

       *subplot*: subplot object
         where to plot SFH

       *legend*: boolean
         whether to draw a legend or not

    By default, sfh will use the formation mass of the star.  In tipsy, this will be
    taken from the starlog file.  Set massform=False if you want the final (observed)
    star formation history

    **Usage:**

    >>> import pynbody.plot as pp
    >>> pp.sfh(s,linestyle='dashed',color='k')


    '''
    if subplot:
        plt = subplot
    else :
        import matplotlib.pyplot as plt

    if 'nbins' in kwargs:
        bins=kwargs['nbins']
        del kwargs['nbins']

    if zrange:
        assert len(zrange) == 2
    else:
        zrange = [0,20]
    if logbins==True:
	logbinsize = (np.log10(zrange[1]+1)-np.log10(zrange[0]+1))/bins
	logzplusonebins = np.log10(zrange[0]+1)+np.arange(bins+1)*logbinsize
	zplusonebins = 10**logzplusonebins
	zbins = zplusonebins-1
    if not logbins:
	zbinsize = (zrange[1]-zrange[0])/bins
	zbins = zrange[0] +np.arange(bins+1)*zbinsize
    from pynbody.analysis import pkdgrav_cosmo as cosmo
    c = cosmo.Cosmology(sim=sim)

    timebins = np.zeros(bins+1)
    for ii in range(bins+1):
	timebins[ii] = 13.7-13.7*c.Exp2Time(1.0 / (1+zbins[ii]))/c.Exp2Time(1)

    print timebins
    binnorm = np.zeros(bins)
    for i in range(bins):
	binnorm[i] = 1e-9 / (timebins[i+1] - timebins[i])

    trangefilt = filt.And(filt.HighPass('tform',str(13.7-timebins[bins])+' Gyr'),
                          filt.LowPass('tform',str(13.7-timebins[0])+' Gyr'))
    tforms = sim.star[trangefilt]['tform'].in_units('Gyr')
    tformslookback = 13.7-tforms

    if massform :
        try:
            weight = sim.star[trangefilt]['massform'].in_units('Msol')# * binnorm / Volume
        except (KeyError, units.UnitsException) :
            warnings.warn("Could not load massform array -- falling back to current stellar masses", RuntimeWarning)
            weight = sim.star[trangefilt]['mass'].in_units('Msol')# * binnorm / Volume
    else:
        weight = sim.star[trangefilt]['mass'].in_units('Msol')# * binnorm / Volume
    if clear : plt.clf()
    sfdens, thebins = np.histogram(tformslookback, weights=weight, bins=timebins)
    sfdens = sfdens*binnorm/Volume
    if not histogram:
	zbincenter = np.zeros(bins)
	for i in range(bins):
		zbincenter[i] = 0.5*(zbins[i]+zbins[i+1])
	print zbincenter
	print sfdens
	print binnorm
	print Volume
	if pltloglog:
		plt.loglog(zbincenter+1,sfdens,label=label,**kwargs)
	else:
		plt.plot(zbincenter,sfdens,label=label,**kwargs)
		plt.yscale('log')
    else:
	sfdens_hist = np.zeros((bins)*2)
        bins_hist = np.zeros((bins+1)*2)
	for i in range(bins*2):
        	sfdens_hist[i] = sfdens[i/2]
        for i in range((bins+1)*2):
                bins_hist[i] = zbins[i/2]
                
        bins_hist = bins_hist[np.arange((bins+1)*2-2)+1]
    	if pltloglog:
		plt.loglog(bins_hist,sfdens_hist,label=label,**kwargs)
	else:
		plt.plot(zbincenter,sfdens,label=label,**kwargs)
                plt.yscale('log')
    if not overplot:
        if not subplot:
                plt.ylim(0.0,1.2*np.max(sfdens))
                plt.xlim(zrange)
                if pltloglog:
			plt.xlabel('log(1+z)',fontsize='large')
		else:
			plt.xlabel('z',fontsize='large')
		plt.ylabel(' Specific SFR [M$_\odot$ yr$^{-1} Mpc^{-3}$]',fontsize='large')
        else:
                plt.set_ylim(0.0,1.2*np.max(sfdens))


def SFR_v_StellarMass(sim,h,dt,grprange=[1,10],grps=False,**kwargs):
	if np.size(grps)<2: nhalos = grprange[1]-grprange[0] + 1
        else: nhalos = np.size(grps)
	SFR = np.array([])
	Mass = np.array([])
	tcurrent = sim.stars['tform'].in_units('Gyr').max()
	if not grps:
		for hh in range(grprange[0],grprange[1]+1):
			if len(h[hh].stars) > .5:
				o, = np.where(tcurrent - h[hh].stars['tform'].in_units("Gyr") < dt)
				massform = array.SimArray(h[hh].stars['massform'][o],sim.stars['mass'].units)
				totalmass = massform.in_units('Msol').sum()
				SFR = np.append(SFR,totalmass/(dt*1e9))
				Mass = np.append(Mass,h[hh].stars['mass'].sum().in_units('Msol'))
	else:
		for hh in range(nhalos):
                        if len(h[grps[hh]].stars) > .5:
                                o = np.where(tcurrent - h[grps[hh]].stars['tform'].in_units("Gyr") < dt)
                                massform = array.SimArray(h[hh].stars['massform'][o],sim.stars['mass'].units)
				totalmass = massform.in_units('Msol').sum()
				SFR = np.append(SFR,totalmass/(dt*1e9))
                                Mass = np.append(Mass,h[grps[hh]].stars['mass'].sum().in_units('Msol'))
	print Mass, SFR
	ind = np.argsort(Mass)
	Mass = Mass[ind]
	SFR = SFR[ind]
	return Mass, SFR

def BHMass_StellarMass(h,grprange=[1,10],grps=False,**kwargs):
	if not grps: nhalos = grprange[1]-grprange[0] + 1
	else: nhalos = np.size(grps)
	Mbh = np.zeros(nhalos)
	Mstar = np.zeros(nhalos)
	cnt = 0
	if not grps:
		for hh in range(grprange[0],grprange[1]+1):
			bhind = np.where(h[hh].stars['tform'] < 0)
			Mbh[cnt] = h[hh].stars['mass'][bhind[0]].in_units('Msol').max()
			Mstar[cnt] = h[hh].stars['mass'].in_units('Msol').sum()
			cnt = cnt + 1
	else:
		for hh in range(nhalos):
			bhind = np.where(h[grp[hh]].stars['tform'] < 0)
                        Mbh[cnt] = h[grp[hh]].stars['mass'][bhind[0]].in_units('Msol').max()
                        Mstar[cnt] = h[grp[hh]].stars['mass'].in_units('Msol').sum()
                        cnt = cnt + 1
			
	ind = np.argsort(Mstar)
	Mstar = Mstar[ind]
	Mbh = Mbh[ind]
	plt.loglog(Mstar,Mbh,**kwargs)	
