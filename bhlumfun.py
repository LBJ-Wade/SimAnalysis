import pynbody
import bhanalysis
import numpy as np
import matplotlib.pyplot as plt
import gc
import readcol
import os
import pickle

z =     np.array([1.  ,1.5 ,2.  ,2.5 ,3.  ,3.5 ,4.  ,4.5 ,5.  ,5.5 ,6.])
zbinsl =np.array([0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75])
zbinsh =np.array([1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25])
c = 2.99792458e10
lbol_sun = 3.9e33
loglbol_sun = np.log10(lbol_sun)
plt.ion()

#Barger_2013 numbers for z = 6 quasars
logphi6B = np.log10(1e-6)
errlogphi6B = 0.6
logLbol6B = 45.05
errlogLbol6B = 0.5

#Fiore+ 2012 number for z>5.8 quasars
logphi6F = np.log10(0.66e-5)
errlogphi6Fp = np.log10(0.66e-5+1.1e-5) - logphi6F
errlogphi6Fm = logphi6F - np.log10(0.66e-5-0.5e-5)
logLbol6F = 45.55
logLbol6Fm = 45.55 - 44.9
logLbol6Fp = 46.2 - 45.55

#convert monochromatic AB absolute magnitude to log(luminosity [ergs/s])
def mcABconv(mag,nu):
        C = 20.638
        return -(2./5.)*mag + C + np.log10(nu)

def mkAbridgeOrbit(simname,s,lmin=1e43,mmin=1e6):
	munits = s.s['mass'].units
	tunits = s.s['x'].units/s.s['vel'].units
	mdotunits = munits/tunits

	Mlimitsim = mmin/munits.in_units('Msol')
        mdotlimit = lmin/(0.1*3e10*3e10)
        mdotlimit /= mdotunits.in_units('g s**-1')
        cstr = """ awk '{if ($4 - $13> """+str(Mlimitsim)+""" && $12 > """+str(mdotlimit)+""") print  $4 " " $12 " " $13 " " $15 " " $16}' """ + simname + ".orbit > " + simname + ".BHorbit.abridged"
	os.system(cstr)
	return

def getLumFun(sim,simname,bins=50,loglmin=43,loglmax=46,vol=25**3,minm=1e6,filename='LumFun.pkl'):
	munits = sim.s['mass'].units
        tunits = sim.s['x'].units/sim.s['vel'].units
        mdotunits = munits/tunits

        tbinsh =np.array([bhanalysis.getTime(zz,sim) for zz in zbinsl])
        tbinsl = np.array([bhanalysis.getTime(zz,sim) for zz in zbinsh])
        dtbins = tbinsh - tbinsl
        dlogl = np.float(loglmax-loglmin)/bins
	if not os.path.exists(simname+'.BHorbit.abridged'):  mkAbridgeOrbit(simname,sim,lmin=10**loglmin,mmin=minm)	
	mass, mdot, dm, dt, scale = readcol.readcol(simname+'.BHorbit.abridged',twod=False)
	ok, = np.where((mass - dm > minm/munits.in_units("Msol"))&(mdot > 10**loglmin/(0.1*3e10*3e10*mdotunits.in_units('g s**-1'))))
	del(dm)
	del(mass)
	gc.collect()
	mdot = pynbody.array.SimArray(mdot[ok],mdotunits)
	dt = pynbody.array.SimArray(dt[ok],tunits)
	scale = scale[ok]
	del(ok)
	gc.collect()
	lum = mdot.in_units('g s**-1')*0.1*3e10*3e10
	del(mdot)
	gc.collect()
	data = np.zeros((len(z),bins))
	for i in range(len(z)):
		print 'redshift ', z[i]
		oo, = np.where((scale**-1 -1 > zbinsl[i])&(scale**-1 -1 < zbinsh[i]))
		weights = dt[oo].in_units('Gyr')/(dtbins[i]*dlogl*vol)
		lumhist, lumbins = np.histogram(np.log10(lum[oo]),range=[loglmin,loglmax],weights=weights,bins=bins)
		data[i,:] = lumhist
		del(oo)
		gc.collect()
	del(lum)
	gc.collect()
	if filename:
		print "saving data..."
		f = open(filename,'wb')
		pickle.dump([data,lumbins],f)
		f.close()
	return data, lumbins

def pltLumFun(data,lumbins,color='blue',linestyle='-',redshift=1,overplot=False,plotdata=True,label=None,linewidth=2):
	zz, = np.where(z==redshift)
	plt.step(lumbins,np.log10(np.append(data[zz,:],data[zz,-1])),color=color,linestyle=linestyle,label=label,lw=linewidth)
	if plotdata==True:
		obs = readcol.readcol('/u/sciteam/tremmel/QSOdata/bol_lf_point_dump.dat',twod=False,asdict=True,skipline=38)
                obs2 = readcol.readcol('/u/sciteam/tremmel/QSOdata/M1450z5_McGreer13.dat',twod=False,asdict=True,skipline=1)
		tt, = np.where(obs['redshift']==redshift)
		plt.errorbar(obs['lbol'][tt] + loglbol_sun, obs['dphi'][tt],yerr=obs['sig'][tt],fmt='o',color='grey',ecolor='grey',label='Hopkins+ 2007 (Compilation)')
		if z[zz] == 6:
			plt.errorbar([logLbol6B],[logphi6B],xerr=errlogLbol6B,yerr=errlogphi6B,fmt='^',color='k',label='Barger+2003')
	                plt.errorbar([logLbol6F],[logphi6F],xerr=[[logLbol6Fm],[logLbol6Fp]],yerr=[[errlogphi6Fm],[errlogphi6Fp]],fmt='s',color='k',label='Fiore+ 2012')
		if z[zz] == 5:
			l1450 = np.log10(4.4)+mcABconv(obs2['M1450'],c/(0.145e-4))
	                dphi = 10**obs2['logphi']
	                dphip = (2./5.) * (dphi+obs2['sig'])
	                dphim = (2./5.) * (dphi - obs2['sig'])
	                dphi = np.log10((2./5.)*dphi)
	                dphierr = [dphi-np.log10(dphim),np.log10(dphip)-dphi]
	                plt.errorbar(l1450,dphi,yerr=dphierr,fmt='D',color='k',label='McGreer+ 2013')
	if overplot==False:
		plt.title(str(zbinsl[zz[0]])+' < z < '+str(zbinsh[zz[0]]))
     		plt.xlabel(r'log$_{10}$($L_{bol}$ [ergs/s]))',fontsize=30)
     		plt.ylabel(r'log$_{10}$($\phi$ [Mpc$^{-3}$ dex$^{-1}$])',fontsize=30)
        plt.legend(loc='lower left',fontsize=20)
	return
