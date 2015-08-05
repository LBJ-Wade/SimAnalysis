import  numpy as np
from math import log
import sys
import pynbody
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import pickle
import os
import gc

def mkDecomp(s,h,minm=10.5,maxm=13,angmom_size="3 kpc"):
        cnt = 1
        while cnt< len(h):
                print "making decomp for halo", cnt
                try:
			Mhalo = h[cnt]['mass'].sum()
		except:
			print "halo", cnt, "does not exist..."
			cnt += 1
			continue
                print np.log10(Mhalo)
                if np.log10(Mhalo) > maxm:
                        print "halo too big"
                        cnt += 1
                        continue
                if np.log10(Mhalo) < minm:
                        print "hit minimum limit on mass!"
                        break
                try:
                        decomp = pynbody.analysis.decomp(h[cnt],angmom_size=angmom_size)
                except:
                        print "decomp failed, moving on"
			cnt += 1
                        continue
		cnt += 1
        print "writing array file..."
        s.write_array('decomp')
        return

def runDecompAll(decomplist,minm=10.5,maxm=13,angmom_size="3 kpc"):
        f = open(decomplist,'r')
        files = f.readlines()
        for i in range(len(files)):
                print "making decomp for ", files[i]
                s = pynbody.load(files[i].strip('\n'))
                h = s.halos()
                s.physical_units()
                mkDecomp(s,h,minm=minm,maxm=maxm,angmom_size=angmom_size)
		del(s)
		del(h)
		gc.collect()
        return

def haloCat(lowz,highz,nhalos=50):
        print "matching halo catalog from ", highz, "to ", lowz
        s1 = pynbody.load(highz)
        s2 = pynbody.load(lowz)
        if s1.properties['a'] > s2.properties['a']:
                print "uh oh! highz file must actually be at higher z!"
                return
        b = pynbody.bridge.OrderBridge(s1,s2)
        cat = b.match_catalog()
        filename = highz+'.cat.z'+str(round(s2.properties['a']**-1-1,3))
        f = open(filename,'wb')
        pickle.dump(cat,f)
        f.close()
	del(s2)
	del(s1)
	gc.collect()
        return

def runhaloCatAll(halocatfile,lowz,nhalos=50):
        r = open(halocatfile,'r')
        files = r.readlines()
        for i in range(len(files)):
                highz = files[i].strip('\n')
                haloCat(lowz,highz,nhalos=nhalos)

        return

def galaxyGrowth(lowz,catfiles,halonum=[1],rmaxlum=None,rinner = 2.5):
	f = open(catfiles,'r')
	files = f.readlines()
	growth = {'Mtot':np.zeros((len(halonum),len(files)+1)),'Mstar':np.zeros((len(halonum),len(files)+1)),'Mgas':np.zeros((len(halonum),len(files)+1)),'MHI':np.zeros((len(halonum),len(files)+1)),'MBH':np.zeros((len(halonum),len(files)+1)),'McenBH':np.zeros((len(halonum),len(files)+1)),'MstarINNER':np.zeros((len(halonum),len(files)+1)),'MgasINNER':np.zeros((len(halonum),len(files)+1)),'MHIINNER':np.zeros((len(halonum),len(files)+1)),'SFR':np.zeros((len(halonum),len(files)+1)),'SFRinner':np.zeros((len(halonum),len(files)+1)),'redshift':np.zeros(len(files)+1)}
	slz = pynbody.load(lowz)
	lz = str(round(slz.properties['a']**-1 -1,3))
	catend = '.cat.z'+lz+'\n'
	for i in range(len(files)):
		xx = files[i].find(catend)
                simname=files[i][0:xx]
		print "getting growth data for", simname
		s = pynbody.load(simname)
                h = s.halos()
		s.physical_units()
		growth['redshift'][i] = s.properties['a']**-1 -1
		catf = open(files[i].strip('\n'))
		cat = pickle.load(catf)
                catf.close()
		for j in range(len(halonum)):
			badcen=0
			print "halo", halonum[j]
                        progs, = np.where(cat==halonum[j])
                        if len(progs)==0:
                                print "no progenitors found in this step!"
                                continue
                        main = progs[0]
			print "progenitor", main
                        h1 = h[main]
			try: pynbody.analysis.halo.center(h1,mode='hyb')
			except: 
				print "cannot find center.. trying to use COM of stars instead..."
				try: pynbody.analysis.halo.center(h1.s,mode='com')
				except:
					print "center failed again"
					badcen = 1
			growth['Mtot'][j,i] = h1['mass'].sum()
			growth['Mstar'][j,i] = h1.s['mass'].sum()
			growth['Mgas'][j,i] = h1.g['mass'].sum()
			growth['MHI'][j,i] = np.sum(h1.g['mass']*h1.g['HI'])
			bhind, = np.where(h1.s['tform']<0)
			if len(bhind)>0:
				growth['MBH'][j,i] = h1.s['mass'][bhind][(h1.s['mass'][bhind]==np.float(h1.s['mass'][bhind].max()))]
				if badcen == 0: growth['McenBH'][j,i] = h1.s['mass'][bhind][(h1.s['r'][bhind]==np.float(h1.s['r'][bhind].min()))]
			growth['SFR'][j,i] = h1.s['massform'][(h1.s['tform'].in_units('Gyr')>(s.properties['time'].in_units('Gyr')-0.1))].in_units('Msol').sum()/1.0e8
			if badcen == 0:
				growth['SFRinner'][j,i] = h1.s['massform'][((h1.s['tform'].in_units('Gyr')>(s.properties['time'].in_units('Gyr')-0.1))&(h1.s['r'].in_units('kpc')<rinner))].in_units('Msol').sum()/1.0e8
				growth['MstarINNER'][j,i] = h1.s['mass'][(h1.s['r'].in_units('kpc')<rinner)].sum()
				growth['MgasINNER'][j,i] = h1.g['mass'][(h1.g['r'].in_units('kpc')<rinner)].sum()
				growth['MHIINNER'][j,i] = np.sum(h1.g['mass'][(h1.g['r'].in_units('kpc')<rinner)]*h1.g['HI'][(h1.g['r'].in_units('kpc')<rinner)])
		del(s)
		del(h)
		del(h1)
		del(cat)
		gc.collect()
	print "calculating values for final step"
        h = slz.halos()
	slz.physical_units()
        for j in range(len(halonum)):
                h1 = h[halonum[j]]
		pynbody.analysis.halo.center(h1,mode='hyb')
                growth['Mtot'][j,i+1] = h1['mass'].sum()
                growth['Mstar'][j,i+1] = h1.s['mass'].sum()
                growth['Mgas'][j,i+1] = h1.g['mass'].sum()
                growth['MHI'][j,i+1] = np.sum(h1.g['mass']*h1.g['HI'])
                bhind, = np.where(h1.s['tform']<0)
		if len(bhind)>0:
                	growth['MBH'][j,i+1] = h1.s['mass'][bhind][(h1.s['mass'][bhind]==np.float(h1.s['mass'][bhind].max()))]
                	growth['McenBH'][j,i+1] = h1.s['mass'][bhind][(h1.s['r'][bhind]==np.float(h1.s['r'][bhind].min()))]
		growth['SFR'][j,i] = h1.s['massform'][(h1.s['tform'].in_units('Gyr')>(slz.properties['time'].in_units('Gyr')-0.1))].in_units('Msol').sum()/1.0e8
		growth['SFRinner'][j,i] = h1.s['massform'][((h1.s['tform'].in_units('Gyr')>(slz.properties['time'].in_units('Gyr')-0.1))&(h1.s['r'].in_units('kpc')<rinner))].in_units('Msol').sum()/1.0e8
                growth['MstarINNER'][j,i+1] = h1.s['mass'][(h1.s['r'].in_units('kpc')<rinner)].sum()
                growth['MgasINNER'][j,i+1] = h1.g['mass'][(h1.g['r'].in_units('kpc')<rinner)].sum()
                growth['MHIINNER'][j,i+1] = np.sum(h1.g['mass'][(h1.g['r'].in_units('kpc')<rinner)]*h1.g['HI'][(h1.g['r'].in_units('kpc')<rinner)])
	return growth

CANDELS_M31 = {'redshift':[0.45,  0.8, 1.0,1.25, 1.55, 1.85,  2.1,  2.5,3.15],
		'lMstar':[10.85,10.81,10.8,10.7,10.62,10.48,10.36,10.15, 9.8],
		'n': [4.2,3.6,3.0,2.5,2.2,1.8,1.0,1.1,1.3],
		'n+':[1.3,1.3,1.2,2.5,1.7,1.7,1.5,2.0,1.9],
		'n-':[1.5,1.0,1.3,1.3,1.3,1.1,0.5,0.6,0.7],
		'UV':[2.0,1.9,1.7,1.7,1.6,1.5,1.2,0.9,0.6],
		'UV+':[0.2,0.2,0.2,0.2,0.3,0.3,0.5,0.6,0.4],
		'UV-':[0.2,0.3,0.3,0.4,0.3,0.3,0.4,0.3,0.3],
		'VJ': [1.3,1.3,1.4,1.2,1.3,1.3,1.1,0.8,0.3],
		'VJ+':[0.1,0.2,0.2,0.3,0.4,0.4,0.5,0.5,0.9],
		'VJ-':[0.1,0.2,0.2,0.2,0.2,0.3,0.4,0.4,0.6]
		}
CANDELS_MW = { 'redshift':[0.45,  0.8,  1.0, 1.25, 1.55,1.85, 2.1, 2.5],
		'lMstar': [10.6,10.47,10.35,10.21,10.06,9.88, 9.7,9.48],
		'n':  [3.4,2.7,2.1,1.5,1.2,1.1,1.3,1.3],
		'n+': [1.7,1.4,1.5,2.1,1.5,1.1,1.4,1.4],
		'n-': [1.7,1.4,1.2,0.8,0.6,0.6,0.6,0.7],
		'UV': [1.9,1.7,1.6,1.4,1.1,0.9,0.6,0.6],
		'UV+':[0.2,0.2,0.3,0.4,0.5,0.4,0.4,0.3],
		'UV-':[0.3,0.3,0.3,0.4,0.3,0.3,0.3,0.3],
		'VJ': [1.3,1.3,1.2,1.2,1.0,0.8,0.5,0.3],
		'VJ+':[0.2,0.3,0.4,0.4,0.4,0.5,0.4,0.5],
		'VJ-':[0.1,0.2,0.2,0.3,0.3,0.4,0.3,0.4]
		}


ABcorr = {'u':0.79,'b':-0.09,'v':0.02,'r':0.21,'i':0.45,'j':0.91,'h':1.39,'k':1.85}
lamcen = {'u':0.365,'b':0.445,'v':0.551,'r':0.658,'i':0.806,'j':1.22,'h':1.63,'k':2.19}

def k(lam,Rv):
	if lam >= 0.63 and lam <= 2.2:
		return 2.659*(1.04/lam - 1.857) + Rv
	if lam >= 0.12 and lam < 0.63:
		return 2.659*(0.011/lam**3 - 0.198/lam**2 + 1.509/lam - 2.156) + Rv
	if lam < 0.12 or lam > 2.2:
		raise ValueError, "wavelength not in acceptable range"

def dustCor(h,s,Rv,A1600max=5.0):
	dustf = 0.01  #Draine 2007
	a = pynbody.array.SimArray(0.1,'1e-6 m') #Todini+Ferrarra (2001), Nozawa+ (2003)
	rho = pynbody.array.SimArray(2.5,'g cm**-3')
#	Rv = 3.1	
	s.physical_units()
	dustExt = {'u':0,'b':0,'v':0,'r':0,'i':0,'j':0,'h':0,'k':0}
	Md = np.sum(h.gas['mass'].in_units('Msol')*(h.g['metals']/0.02)*h.gas['HI']*dustf)
#	Md = np.sum(h.gas['mass'].in_units('Msol')*(h.g['OxMassFrac']/(h.g['hydrogen']*16.*2.5e-4))*h.gas['HI']*dustf)
	#Md = np.sum(h.gas['mass'].in_units('Msol')*(10**h.g['oxh'])*h.gas['HI']*dustf)
	try:
		pynbody.analysis.halo.center(h,vel=True,mode='hyb',wrap=True)
	except:
		 pynbody.analysis.halo.center(h,vel=False,mode='hyb',wrap=True)
#	Rhl = pynbody.analysis.luminosity.half_light_r(h.s[(h.s['tform']>0)],band='v')
	#Rgas = h.gas['r'][(h.gas['HI']>0.5)].max()
#	rhalo = h.dm['r'].in_units('kpc').max()
	rord = np.argsort(h.gas['r'])
	Msum = np.cumsum(h.gas['mass'][rord].in_units('Msol')*(h.g['metals'][rord]/0.02)*h.gas['HI'][rord]*dustf)
#	Msum = np.cumsum(h.gas['mass'][rord].in_units('Msol')*(10**h.g['oxh'][rord])*h.gas['HI'][rord]*dustf)
#	Msum = np.cumsum(h.gas['mass'][rord].in_units('Msol')*(h.g['OxMassFrac'][rord]/(h.g['hydrogen'][rord]*16.*2.5e-4))*h.gas['HI'][rord]*dustf)
	xx, = np.where(Msum >= Md*0.5)
	Rhalf = h.gas['r'][rord[xx]].in_units('kpc')[0]
	#h.g['dustM'] = h.g['mass']*h.g['HI']*(h.g['metals']/0.02)*dustf
	#p = pynbody.analysis.profile.Profile(h.g,nbins=100,type='log')
	#dr = np.append(p['rbins'][1:].in_units('kpc'),h.g['r'].in_units('kpc').max()) - p['rbins'][0:].in_units('kpc')
	#o, = np.where(p['mass'] > 0)
	#sigD = np.sum(p['mass'][o]*p['oxh'][o]*p['HI'][o]*dustf/(4.*np.pi*p['rbins'][o]**2))
	#sigD = np.sum(p['dustM'][o]/(4.*np.pi*p['rbins'][o].in_units('kpc')**2))
	sigD = 0.5*Md/(np.pi*Rhalf**2)
	tau = 3.*sigD/(4.*a.in_units('kpc')*rho.in_units('Msol kpc**-3'))
	A1600 = 1.086*tau
	if A1600 > A1600max: A1600=A1600max
	EBV = A1600/k(0.16,Rv)

	print tau,A1600,Rhalf
	
	for key in dustExt.keys():
		dustExt[key] = k(lamcen[key],Rv)*EBV

	return dustExt,Rhalf

def getAllDust(lowz,catfiles,Rv,A1600max=2.0,halonum=[1],filename='dust.pkl'):
	f = open(catfiles,'r')
	files = f.readlines()
	dustExt = {'halos':halonum,'Rhl':np.zeros((len(halonum),len(files)+1)),'RDhalf':np.zeros((len(halonum),len(files)+1)),'z':np.zeros(len(files)+1),'u':np.zeros((len(halonum),len(files)+1)),'b':np.zeros((len(halonum),len(files)+1)),'v':np.zeros((len(halonum),len(files)+1)),'r':np.zeros((len(halonum),len(files)+1)),'i':np.zeros((len(halonum),len(files)+1)),'j':np.zeros((len(halonum),len(files)+1)),'h':np.zeros((len(halonum),len(files)+1)),'k':np.zeros((len(halonum),len(files)+1))}
	slz = pynbody.load(lowz)
        lz = str(round(slz.properties['a']**-1 -1,3))
        catend = '.cat.z'+lz+'\n'
        for i in range(len(files)):
                print "calculating dust corrections for ", files[i].strip('\n')
                xx = files[i].find(catend)
                simname=files[i][0:xx]
                s = pynbody.load(simname)
                h = s.halos()
                dustExt['z'][i] = s.properties['a']**-1 -1
                s.physical_units()
                catf = open(files[i].strip('\n'))
                cat = pickle.load(catf)
                catf.close()
		for j in range(len(halonum)):
                        print "halo", halonum[j]
                        progs, = np.where(cat==halonum[j])
                        if len(progs)==0:
                                print "no progenitors found in this step!"
                                continue
                        main = progs[0]
                        print "progenitor", main
                        h1 = h[main]
			dust,Rhalf = dustCor(h1,s,Rv,A1600max=A1600max)
			dustExt['RDhalf'][j,i] = Rhalf
#			dustExt['Rhl'][j,i] = Rhl
			for key in dust.keys():
				dustExt[key][j,i] = dust[key]
		del(s)
                del(h)
                del(h1)
                del(cat)
                gc.collect()
	print "calculating values for final step"
        h = slz.halos()
        slz.physical_units()
        for j in range(len(halonum)):
                h1 = h[halonum[j]]
		dust,Rhalf= dustCor(h1,slz,Rv)
		dustExt['z'][i+1] = slz.properties['a']**-1 -1
		dustExt['RDhalf'][j,i+1] = Rhalf
		for key in dust.keys():
			dustExt[key][j,i+1] = dust[key]
	if filename:
                print "saving data..."
                f = open(filename,'wb')
                pickle.dump(dustExt,f)
                f.close()
	return dustExt
	
def colorHistory(lowz,catfiles,halonum=[1],maxr = None,b_band='u',g_band='v',r_band='j',filename='mags.pkl'):
	f = open(catfiles,'r')
	files = f.readlines()
	magnitudes = {'halos':halonum,b_band:np.zeros((len(halonum),len(files)+1)),g_band:np.zeros((len(halonum),len(files)+1)),r_band:np.zeros((len(halonum),len(files)+1)),'z':np.zeros(len(files)+1)}
#	color1 = np.zeros((len(halonum),len(files)+1))
#	color2 = np.zeros((len(halonum),len(files)+1))
#	redshift = np.zeros(len(files))
	slz = pynbody.load(lowz)
	lz = str(round(slz.properties['a']**-1 -1,3))
	catend = '.cat.z'+lz+'\n'
	for i in range(len(files)):
		print "calculating magnitudes for ", files[i].strip('\n')
		xx = files[i].find(catend)
		simname=files[i][0:xx]
		s = pynbody.load(simname)
		h = s.halos()
		magnitudes['z'][i] = s.properties['a']**-1 -1
		s.physical_units()
		catf = open(files[i].strip('\n'))
		cat = pickle.load(catf)
		catf.close()
		for j in range(len(halonum)):
			print "halo", halonum[j]
			progs, = np.where(cat==halonum[j])
			if len(progs)==0:
				print "no progenitors found in this step!"
				continue
			main = progs[0]
			print "progenitor", main
			h1 = h[main]
			if maxr:
				pynbody.analysis.halo.center(h1.s,mode='com')
				use, = np.where((h1.s['r']<maxr)&(h1.s['tform']>0))
			else:
				use, = np.where(h1.s['tform']>0)
			magnitudes[b_band][j,i] = pynbody.analysis.luminosity.halo_mag(h1.s[use],band=b_band) + ABcorr[b_band]
			magnitudes[g_band][j,i] = pynbody.analysis.luminosity.halo_mag(h1.s[use],band=g_band) + ABcorr[g_band]
			magnitudes[r_band][j,i] = pynbody.analysis.luminosity.halo_mag(h1.s[use],band=r_band) + ABcorr[r_band]
			#color1[j,i] = bmag-gmag
			#color2[j,i] = gmag-rmag
		del(s)
		del(h)
		del(h1)
		del(cat)
		gc.collect()
	print "calculating values for final step"
	h = slz.halos()
	slz.physical_units()
	magnitudes['z'][i+1] = slz.properties['a']**-1 -1
	for j in range(len(halonum)):
		h1 = h[halonum[j]]
		if maxr:
                        pynbody.analysis.halo.center(h1.s,mode='com')
                        use, = np.where((h1.s['r']<maxr)&(h1.s['tform']>0))
                else:
                        use, = np.where(h1.s['tform']>0)
                magnitudes[b_band][j,i+1] = pynbody.analysis.luminosity.halo_mag(h1.s[use],band=b_band) + ABcorr[b_band]
                magnitudes[g_band][j,i+1] = pynbody.analysis.luminosity.halo_mag(h1.s[use],band=g_band) + ABcorr[g_band]
                magnitudes[r_band][j,i+1] = pynbody.analysis.luminosity.halo_mag(h1.s[use],band=r_band) + ABcorr[r_band]
		#	color1[j,i+1] = bmag-gmag
		#	color2[j,i+1] = gmag-rmag
	if filename:
		print "saving data..."
		f = open(filename,'wb')
		pickle.dump(magnitudes,f)
		f.close()
	return magnitudes	
