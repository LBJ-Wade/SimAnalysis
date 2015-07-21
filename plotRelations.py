import  numpy as np
from math import log
import sys
import pynbody
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import pickle
import os
import scipy
plt.ion()
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
plt.rc('font', weight='medium')
plt.rc('axes', linewidth=2)
plt.rc('xtick.major',width=2)
plt.rc('ytick.major',width=2)

def moster13allvar(logM,z,M10,M11,N10,N11,b10,b11,g10,g11):
	r = z/(z+1)
	logM1 = M10 + M11 * r
	N = N10 + N11 * r
	b = b10 + b11 * r
	g = g10 + g11 * r
	x = logM - logM1
	ratio = 2 * N / ( 10**(-b * x) + 10**(g * x) )
	return ratio

def moster13(logM, z):
    M10 = 11.590
    M11 =  1.195
    N10 =  0.0351
    N11 = -0.0247
    b10 =  1.376
    b11 = -0.826
    g10 =  0.608
    g11 =  0.329
    return moster13allvar(logM,z,M10,M11,N10,N11,b10,b11,g10,g11)

def partial_derivative(func, var=0, point=[]):
    args = point[:]
    def wraps(x):
        args[var] = x
        return func(*args)
    return scipy.misc.derivative(wraps, point[var], dx = 1e-8)

def errmoster13(logM,z):
	M10 = 11.590
	M11 =  1.195
	N10 =  0.0351
	N11 = -0.0247
	b10 =  1.376
	b11 = -0.826
	g10 =  0.608
	g11 =  0.329
	sigM10 = 0.236
	sigM11 = 0.353
	sigN10 = 0.0058
	sigN11 = 0.0069
	sigb10 = 0.153
	sigb11 = 0.225
	sigg10 = 0.059
	sigg11 = 0.173
	sigvar = [sigM10,sigM11,sigN10,sigN11,sigb10,sigb11,sigg10,sigg11]
	sigma = np.zeros(len(logM))
	for i in range(len(logM)):
		point = [logM[i],z,M10,M11,N10,N11,b10,b11,g10,g11]
		for j in range(8):
			sigma[i] += partial_derivative(moster13allvar,var=j+2,point=point)**2 * sigvar[j]**2
	sigma = np.sqrt(sigma)
	return sigma


def behroozi13(logM,a):
	z = a**-1 - 1
	v = np.exp(-4.*a**2)
	le = -1.777 + (-0.006*(a-1)*v) - 0.119*(a-1)
	lM1 = 11.514 + (-1.793*(a-1)-0.251*z)*v
	A = -1.412+0.731*(a-1)*v
	d = 3.508+(2.608*(a-1)-0.043*z)*v
	g = 0.316+(1.319*(a-1)+0.279*z)*v
	def f(x):
		return -1.0*np.log10(10**(A*x)+1) + d*((np.log10(1+np.exp(x)))**g)/(1+np.exp(10**(-x)))
	lMstar = le +lM1 + f(logM-lM1) - f(0)
	ratio = 10**(lMstar - logM)
	return ratio

def SpecAngMom(beta,logMstar):
        '''
        calculate the predicted angular momentum from Mstar and B/T ratio
        From Obreschkow and Glazebrook, 2013
        '''
        k = 0.89
        a = 0.94
        g = 7.03
        predlogj = log(k, 10) - g * beta + a * (logMstar - 10) + 3
        return predlogj
def beta(logjOverM):
        '''
        calculate Beta from spec ang mom per unit mass.
        From Obreschkow and Glazebrook, 2013
        '''
        k1 = -.3
        k2 = -.01
        a = 1e-7
        beta  = k1*(logjOverM-log(a,10)) + k2
        return beta

def BHMstar(logMstar):
        '''
        predict BH mass given a stellar mass
        based on analysis by Haring and Rix 2004 and Schramm + Silverman 2013
        '''
        c = 8.31
        a = 1.12
        b = 11.
        predlogMBH = c + a * (logMstar-b)
        return predlogMBH
def BHMBulge(logMbulge):
        '''
        predict BH mass given a bulge mass
        based on analysis by Haring and Rix 2004 and Kormendy and Ho 2013
        '''
        c = 8.69
        a = 1.16
        b = 11.
        predlogMBH = c + a * (logMbulge-b)
        return predlogMBH


def HIFrac(logMstar):
        '''
        predict HI mass fraction given a stellar mass
        based on data from SHIELD and ALFALFA
        '''
        a = 5.0408
        b = 0.5404
        predlogFg = a - b * logMstar
        return predlogFg


def plotSMHM(s,h,minm=10,maxm=13,skiphalo=[],plottype='Ratio',plotfit=['mos','beh'],plotfitz=[0,0.5],correct=True,findSats=True,Satlist=[],pntstyle='go',satpntstyle='b*',msize=20,lnstyle=[['b-','r-'],['b--','r--']],lnthick=[1.5,1.5],axes = None,label=None,legend=True,filename='SMHM.pkl'):
	'''
	plot SMHM for simulation s for halos in h with total (log) virial mass in range[minm,maxm]
	plottype = Ratio or Mass if you want the y axis to be Mstar/Mtot or Mstar
	plotfit = True/False to plotting the observed relation (Moster 13)
	correct = True/False to applying the Munshi 13 correction to stellar mass and halo mass (Mstar x 0.6 and Mhalo/0.8)
	pntstype = the stype of the data points on the plot. Same syntax as matplotlib, e.g. 'go' means green circles, 'bx' means blue x's
	lncolor = the color of the fit line
	axes = the axes object you want the plot drawn in
	'''
	s.physical_units()
	cnt = 1
	lMhalo = np.array([])
	lMstar = np.array([])
	sats = np.array([])
	if os.path.exists(filename):
		f = open(filename)
		lMhalo,lMstar,sats = pickle.load(f)
		f.close()
	else:
		while cnt < len(h):
			print "getting data for halo ", cnt
			if len(skiphalo)>0: 
				bad, = np.where(skiphalo==cnt)
				cnt += 1
				if len(bad) > 0:
					print "skipping this one!" 
					continue
			Mhalo = h[cnt]['mass'].sum()
			if np.log10(Mhalo) > maxm: 
				cnt += 1
				continue
	                if np.log10(Mhalo) < minm: 
				print "hit minimum limit on mass!"
				break
			Mstar = h[cnt].stars['mass'].sum()
			if findSats == True and not Satlist:
				print "determining if satellite..."
				pynbody.analysis.halo.center(h[cnt],wrap=True,mode='hyb',vel=False)	
				r = h[cnt]['r'].max()
				grps = s['amiga.grp'][((s['r']>r)&(s['r']<1.5*r))]
				ugrps,count = np.unique(grps,return_counts=True)
				oo, = np.where(count.astype(np.float)/len(grps)>0.5)
				satflag = 0
				if len(oo)>0:
					if (ugrps[oo] < cnt) and (ugrps[oo] > 0):
						print "Yes"
						satflag=1
					else: 
						satflag=0
						print "no"
				if satflag == 0:
					if count[((ugrps<cnt)&(ugrps>0))].astype(np.float).sum()/len(grps) >0.7:
						satflag=1
						print "yes"
					else:
						satflag=0
					
				sats = np.append(sats,satflag)	
			if len(Satlist)>0:
				oss, = np.where(np.array(Satlist)==cnt)
				if len(oss)>0: sats = np.append(sats,1)
				else: sats = np.append(sats,0)
			if Mstar == 0: 
				print "halo has no stars... skipping"
				cnt += 1
				continue
			if correct==True: Mhalo /= 0.8
			lMhalo = np.append(lMhalo,np.log10(Mhalo))
			if correct==True: Mstar *= 0.6
			lMstar = np.append(lMstar,np.log10(Mstar))
			cnt += 1
		print sats
	if plotfit:
		lmhaloline = np.arange(minm-1.0,maxm+0.5,0.01)
		ratios = {'label':[]}
		sigma = {}
		ratios['label'] = np.zeros((len(plotfit),len(plotfitz))).astype(np.str)
		for fit in plotfit:
			sigma[fit] = np.zeros((len(plotfitz),len(lmhaloline)))
			ratios[fit]=np.zeros((len(plotfitz),len(lmhaloline)))
		print ratios
		for i in range(len(plotfit)):
			for zz in range(len(plotfitz)):
 		       		if plotfit[i] == 'mos': 
					ratios[plotfit[i]][zz,:] = moster13(lmhaloline,plotfitz[zz])
					ratios['label'][i,zz] = "Moster+ 13, z ="+str(plotfitz[zz])
					sigma[plotfit[i]][zz,:] = errmoster13(lmhaloline,plotfitz[zz])
				if plotfit[i] == 'beh': 
					ratios[plotfit[i]][zz,:] = behroozi13(lmhaloline,plotfitz[zz])
					ratios['label'][i,zz]  = "Behroozi+ 13, z ="+str(plotfitz[zz])
	if findSats == False and len(Satlist)==0:
		sats = np.zeros(len(lMstar))
	if plottype == 'Ratio':
		if not axes:
			if satpntstyle:
				plt.plot(lMhalo[(sats==0)],lMstar[(sats==0)]-lMhalo[(sats==0)],pntstyle,markersize=msize,label=label+' (central)')
				plt.plot(lMhalo[(sats==1)],lMstar[(sats==1)]-lMhalo[(sats==1)],satpntstyle,markersize=msize,label=label+' (sat/int)')
			else:
				plt.plot(lMhalo,lMstar-lMhalo,pntstyle,label=label,markersize=msize)
			if plotfit:
				for i in range(len(plotfit)):
					for j in range(len(plotfitz)):
						plt.plot(lmhaloline,np.log10(ratios[plotfit[i]][j,:]),lnstyle[i][j],label=ratios['label'][i,j],linewidth=lnthick[i])
						plt.fill_between(lmhaloline,np.log10(ratios[plotfit[i]][j,:]-sigma[plotfit[i]][j,:]),np.log10(ratios[plotfit[i]][j,:]+sigma[plotfit[i]][j,:]),facecolor='grey',alpha=0.5)
				
		else:
                        axes.plot(lMhalo[(sats==0)],lMstar[(sats==0)]-lMhalo[(sats==0)],pntstyle,label=label+' (central)')
                        axes.plot(lMhalo[(sats==1)],lMstar[(sats==1)]-lMhalo[(sats==1)],satpntstyle,label=label+' (sat/int)')
			if plotfit:
				for i in range(len(plotfit)):
					axes.plot(lmhaloline,np.log10(ratios[plotfit[i]]),lnstyle[i][j],label=ratios['label'][i],linewidth=lnthick[i])

	if plottype == 'Mass':
		if not axes: 
                        plt.plot(lMhalo,lMstar,pntstyle,label=label)
			if plotfit:
				for i in range(len(plotift)):
                        	        plt.plot(lmhaloline,np.log10(ratios[plotfit[i]]*10**lmhaloline),lnstyle[i],label=ratios['label'][i],linewidth=lnthick[i])
                else:   
                        axes.plot(lMhalo,lMstar-lMhalo,pntstyle,label=label)
			if plotfit:
				for i in range(len(plotift)):
                        		axes.plot(lmhaloline,np.log10(ratios[plotfit[i]]*10**lmhaloline),lnstyle[i],label=ratios['label'][i],linewidth=lnthick[i])
	if legend:
		if not axes: plt.legend()
		else: axes.legend()
	if filename:
		f = open(filename,'wb')
		pickle.dump([lMhalo,lMstar,sats],f)
		f.close()
	return

def mkDecomp(s,h,minm=10.5,maxm=13,angmom_size="3 kpc"):
	cnt = 1
	while cnt< len(h):
		print "making decomp for halo", cnt
		Mhalo = h[cnt]['mass'].sum()
		print np.log10(Mhalo)
		if np.log10(Mhalo) > maxm:
			print "halo too big"
                        cnt += 1
                        continue
                if np.log10(Mhalo) < minm:
                        print "hit minimum limit on mass!"
                        break 
		
		pynbody.analysis.halo.center(h[cnt],mode='hyb',wrap=True)
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
	return

def runhaloCatAll(halocatfile,lowz,nhalos=50):
	r = open(halocatfile,'r')
	files = r.readlines()
	for i in reante(len(files)):
		highz = files[i].strip('\n')
		haloCat(lowz,highz,nhalos=nhalos)
		
	return

#def isSat(s,h,rmax=2.0,fracth=0.5):
#	rad = h[r].in_units('kpc').max()
#	near = s[(s

def plotBHStar(s,h,minm=10.5,maxm=13,skiphalo=[],BHtype='mass', rmax=2, plotfit=True,pntstyle='go',lnstyle='k-',msize=10,axes = None,label=None,legend=True,filename='BHMstar.pkl'):
	s.physical_units()
        cnt = 1
	lMStar = np.array([])
        lMBH = np.array([])
	if filename and os.path.exists(filename):
		f = open(filename,'rb')
		lMStar,lMBH = pickle.load(f)
		f.close()
	if not filename or not os.path.exists(filename):
		while cnt < len(h):
	                print "getting data for halo ", cnt
	                if len(skiphalo)>0:
	                        bad, = np.where(skiphalo==cnt)
	                        cnt += 1
	                        if len(bad) > 0:
	                                print "skipping this one!"
	                                continue
	                Mhalo = h[cnt]['mass'].sum()
	                if np.log10(Mhalo) > maxm:
	                        cnt += 1
	                        continue
	                if np.log10(Mhalo) < minm:
	                        print "hit minimum limit on mass!"
	                        break
	                Mstar = h[cnt].stars['mass'].sum()
	                if Mstar == 0:
	                        print "halo has no stars... skipping"
	                        cnt += 1
	                        continue
			pynbody.analysis.halo.center(h[cnt],mode='hyb',wrap=True)
	                bhs, = np.where((h[cnt].stars['tform']<0)&(h[cnt].stars['r']<rmax))
        	        if len(bhs)==0:
        	                print "No BH in this halo"
        	                cnt +=1
               	        	continue
                	if h[cnt].stars['r'][bhs].min() > rmax:
                	        print "No BH near the center of this halo"
                	        cnt += 1
                	        continue
                	if BHtype=='mass':
                	        Mbh = h[cnt].stars['mass'][bhs][(h[cnt].stars['mass'][bhs]==np.float(h[cnt].stars['mass'][bhs].max()))]
                	if BHtype=='center':
                	        Mbh = h[cnt].stars['mass'][bhs][(h[cnt].stars['r'][bhs]==np.float(h[cnt].stars['r'][bhs].min()))]
                	if BHtype=='bulge':
				bulge = h[cnt].s[h[cnt].s['decomp'] == 3]
                	        rmaxb = bulge['r'].max()
                	        if h[cnt].stars['r'][bhs].min() > rmaxb:
                	                print "halo has no BHs with Bulge region"
                	                continue
                	        Mbh = h[cnt].stars['mass'][bhs][(h[cnt].stars['r'][bhs]<rmaxb)].max()
			lMBH = np.append(lMBH,np.log10(Mbh))
                	lMStar = np.append(lMStar,np.log10(Mstar))
                	cnt += 1
	if filename and not os.path.exists(filename):
		f = open(filename,'wb')
		pickle.dump([lMStar,lMBH],f)
		f.close()
        if plotfit:
                lmstarline = np.arange(9,lMStar.max()+1,0.1)
                bhmassline = BHMstar(lmstarline)
        if not axes:
                        plt.plot(lMStar,lMBH,pntstyle,label=label,markersize=msize)
                        if plotfit: 
				plt.plot(lmstarline,bhmassline,lnstyle,label='Schramm + Silverman 2013')
				plt.fill_between(lmstarline,bhmassline-0.3,bhmassline+0.3,facecolor='grey',alpha = 0.25)
        else:
                        axes.plot(lMStar,lMBH,pntstyle,label=label,markersize=msize)
                        if plotfit: 
				axes.plot(lmstarline,bhmassline,lnstyle,lnstyle,label='Schramm + Silverman 2013')
				axes.fill_between(lmstarline,bhmassline-0.3,bhmassline+0.3,facecolor='grey',alpha = 0.25)

        if legend:
                if not axes: plt.legend()
                else: axes.legend()
        if not plotfit: return lMStar, lMBH
	if plotfit: return lMStar, lMBH,lmstarline,bhmassline
def plotBHBulge(s,h,minm=10.5,maxm=13,skiphalo=[],BHtype='mass', rmax=2, plotfit=True,pntstyle='go',msize=10,lnstyle='k-',axes = None,label=None,legend=True,filename='BHMBulge.pkl'):
	'''
        plot SMHM for simulation s for halos in h with total (log) virial mass in range[minm,maxm]
        plotfit = True/False to plotting the observed relation (Moster 13)
        pntstype = the stype of the data points on the plot. Same syntax as matplotlib, e.g. 'go' means green circles, 'bx' means blue x's
        lncolor = the color of the fit line
        axes = the axes object you want the plot drawn in
	BHtype = 'mass', 'central', or 'bulge' which means the code chooses the most massive, most central, or most massive within the extend of the bulge for the analysis
	rmax = the maximum radius allowed for a BH to be counted as "within" the galaxy. in kpc
        '''
        s.physical_units()
        cnt = 1
        lMbulge = np.array([])
        lMBH = np.array([])
	if filename and os.path.exists(filename):
		f = open(filename)
		lMbulge,lMBH = pickle.load(f)
		f.close()
	if not filename or not os.path.exists(filename):
	        while cnt < len(h):
	                print "getting data for halo ", cnt
	                if len(skiphalo)>0:
	                        bad, = np.where(skiphalo==cnt)
	                        cnt += 1
	                        if len(bad) > 0:
					print "skipping this one!" 
					continue
	                Mhalo = h[cnt]['mass'].sum()
	                if np.log10(Mhalo) > maxm:
	                        cnt += 1
	                        continue
	                if np.log10(Mhalo) < minm:
	                        print "hit minimum limit on mass!"
	                        break
	                Mstar = h[cnt].stars['mass'].sum()
	                if Mstar == 0:
	                        print "halo has no stars... skipping"
	                        cnt += 1
	                        continue
	                disk_thin    = h[cnt].s[h[cnt].s['decomp'] == 1]
	                disk_thick   = h[cnt].s[h[cnt].s['decomp'] == 4]
	                bulge        = h[cnt].s[h[cnt].s['decomp'] == 3]
	                bulge_pseudo = h[cnt].s[h[cnt].s['decomp'] == 5]
	                st           = h[cnt].s[h[cnt].s['decomp'] != 2]
	
			if bulge['mass'].sum() == 0:
				print "No Bulge"
				cnt +=1
				continue
			pynbody.analysis.halo.center(h[cnt],mode='hyb',wrap=True)
			bhs, = np.where((h[cnt].stars['tform']<0)&(h[cnt].stars['r']<rmax))
			if len(bhs)==0:
				print "No BH in this halo"
				cnt +=1
				continue
			if h[cnt].stars['r'][bhs].min() > rmax:
				print "No BH near the center of this halo"
				cnt += 1
				continue
			if BHtype=='mass': 
				Mbh = h[cnt].stars['mass'][bhs][(h[cnt].stars['mass'][bhs]==np.float(h[cnt].stars['mass'][bhs].max()))]
			if BHtype=='center':
				Mbh = h[cnt].stars['mass'][bhs][(h[cnt].stars['r'][bhs]==np.float(h[cnt].stars['r'][bhs].min()))]
			if BHtype=='bulge':
				rmaxb = bulge['r'].max()
				if h[cnt].stars['r'][bhs].min() > rmaxb:
					print "halo has no BHs with Bulge region"
					continue
				Mbh = h[cnt].stars['mass'][bhs][(h[cnt].stars['r'][bhs]<rmaxb)].max()
	                lMBH = np.append(lMBH,np.log10(Mbh))
			lMbulge = np.append(lMbulge,np.log10(bulge['mass'].sum()))
	                cnt += 1
	if filename and not os.path.exists(filename):
		f = open(filename,'wb')
		pickle.dump([lMbulge,lMBH],f)
		f.close()
	if plotfit==True:
		lmbulgeline = np.arange(7.5,lMbulge.max()+1,0.1)
		bhmassline = BHMBulge(lmbulgeline)
	if not axes:
                        plt.plot(lMbulge,lMBH,pntstyle,label=label,markersize=msize)
                        if plotfit==True: 
				plt.plot(lmbulgeline,bhmassline,lnstyle,label='Kormendy + Ho 2013')
				plt.fill_between(lmbulgeline,bhmassline-0.3,bhmassline+0.3,facecolor='grey',alpha = 0.25)
        else:
                        axes.plot(lMbulge,lMBH,pntstyle,label=label,markersize=msize)
                        if plotfit==True: 
				axes.plot(lmbulgeline,bhmassline,lnstyle,lnstyle,label='Kormendy + Ho 2013')
				axes.fill_between(lmbulgeline,bhmassline-0.3,bhmassline+0.3,facecolor='grey',alpha = 0.25)

        if legend:
                if not axes: plt.legend()
                else: axes.legend()
        if plotfit==True: return lMbulge, lMBH,lmbulgeline,bhmassline
	else: return lMbulge, lMBH
