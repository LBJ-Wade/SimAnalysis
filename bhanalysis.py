"""

readBHdata
=====

by MJT
"""
import numpy as np
import matplotlib.pyplot as plt
import pynbody
from pynbody.analysis import profile, angmom, halo
from pynbody import snapshot, filt, units, config, array
from scipy import optimize as opt
import warnings
import math
import os
from pynbody.analysis import pkdgrav_cosmo as cosmo
import pickle
import readcol
import gc

plt.ion()
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
plt.rc('font', weight='medium')
plt.rc('axes', linewidth=2)
plt.rc('xtick.major',width=2)
plt.rc('ytick.major',width=2)


def writeBHMark(simname,step,Name=None,iord=False,massrange=False):
	if not Name: f = open('BH.'+step+'.mark','w')
	else: f = open(Name,'w')
	s = pynbody.load(simname+'.'+step)
	f.write(str(len(s))+' '+str(len(s.gas))+' '+str(len(s.star))+'\n')
	if not iord: 
		if not massrange:
			bhind, = np.where(s.stars['tform']<0)
		else:
			if len(massrange) != 2:
				print "error massrange must be a length 2 tuple!"
				return
			bhind, = np.where((s.stars['tform']<0)&(s.stars['mass'].in_units('Msol')<massrange[1])&(s.stars['mass'].in_units('Msol')>massrange[0]))
	else:	
		bhind = np.array([])
		for ii in range(len(iord)):
			tmpind, = np.where(s.stars['iord']==iord[ii])
			if len(tmpind)==0: print "uh oh... iord ", iord[ii]," not found!"
			bhind = np.append(bhind,tmpind)
	bhindreal = bhind+len(s.dark)+len(s.gas)+1
	for ii in range(len(bhindreal)):
		f.write(str(bhindreal[ii])+'\n')
	f.close()
	del(s)
	return bhindreal

def getTime(z,sim):
	c = cosmo.Cosmology(sim=sim)
	return 13.7*c.Exp2Time(1.0 / (1+z))/c.Exp2Time(1)

def getFileLists(simname):
	simname_split = simname.split('.')
	num = len(simname_split)
	os.system('ls  *.iord | cut -d"." -f1-'+str(num+1)+' > files.list')
        os.system('ls *.iord | cut -d"." -f'+str(num+1)+' > steps.list')
	
def trackIsoBH(simname,SF=False,filename=False,BeLazy=False):
        if not os.path.exists('files.list'):
                getFileLists(simname)
	f = open('files.list')
	files = f.readlines()
	distp = array.SimArray(np.zeros(len(files)),'pc')
	xp = array.SimArray(np.zeros(len(files)),'pc')
	yp = array.SimArray(np.zeros(len(files)),'pc')
	zp = array.SimArray(np.zeros(len(files)),'pc')
	vx = array.SimArray(np.zeros(len(files)),'km s**-1')
	vy = array.SimArray(np.zeros(len(files)),'km s**-1')
	vz = array.SimArray(np.zeros(len(files)),'km s**-1')
	vcx = array.SimArray(np.zeros(len(files)),'km s**-1')
        vcy = array.SimArray(np.zeros(len(files)),'km s**-1')
        vcz = array.SimArray(np.zeros(len(files)),'km s**-1')
	vrp = array.SimArray(np.zeros(len(files)),'km s**-1')
	distc = array.SimArray(np.zeros(len(files)),'pc')
        xc = array.SimArray(np.zeros(len(files)),'pc')
        yc = array.SimArray(np.zeros(len(files)),'pc')
        zc = array.SimArray(np.zeros(len(files)),'pc')
	vrc = array.SimArray(np.zeros(len(files)),'km s**-1')
	cnt = 0
	time = array.SimArray(np.zeros(len(files)),'Gyr')
#	time = (np.arange(len(files))+1)*outInterval*dt
        for sim in files:
		print "getting data from ", sim.strip('\n')
        	s = pynbody.load(sim.strip('\n'))
		if BeLazy==False: 
			cen = pynbody.analysis.halo.shrink_sphere_center(s)
			vcen = pynbody.analysis.halo.center_of_mass_velocity(s) 
		cenpot = s['pos'][(s['phi']==float(s['phi'].min()))]
		cenpot = cenpot[0]
		time[cnt] = s.properties['time'].in_units('Gyr')
		if SF==False:
			if BeLazy==False:
				vx[cnt] = s.stars['vx'].in_units('km s**-1')[0]
				vy[cnt] = s.stars['vx'].in_units('km s**-1')[0]
				vz[cnt] = s.stars['vx'].in_units('km s**-1')[0]
				s.stars['vel'] -= vcen
				vcx[cnt] = s.stars['vx'].in_units('km s**-1')[0]
                        	vcy[cnt] = s.stars['vy'].in_units('km s**-1')[0]
                        	vcz[cnt] = s.stars['vz'].in_units('km s**-1')[0] 
			s.stars['pos'] -= cenpot
			xp[cnt] = s.stars['x'].in_units('pc')[0]
			yp[cnt] = s.stars['y'].in_units('pc')[0]
			zp[cnt] = s.stars['z'].in_units('pc')[0]
			if BeLazy==False: vrp[cnt] = s.stars['vr'].in_units('km s**-1')[0]
			distp[cnt] = s.stars['r'].in_units('pc')[0]
                        s.stars['pos'] += cenpot
			if BeLazy==False:s.stars['pos'] -= cen
			xc[cnt] = s.stars['x'].in_units('pc')[0]
                        yc[cnt] = s.stars['y'].in_units('pc')[0]
                        zc[cnt] = s.stars['z'].in_units('pc')[0]
			if BeLazy==False: 
				vrc[cnt] = s.stars['vr'].in_units('km s**-1')[0]	
                        	distc[cnt] = s.stars['r'].in_units('pc')[0]
		cnt += 1
		del(s)
		gc.collect()
	
	BHorbitInfo = {'xp':xp, 'yp':yp, 'zp':zp, 'xc':xc, 'yc':yc, 'zc':zc, 'vrp':vrp, 'vrc':vrc, 'vx':vx, 'vy':vy, 'vz':vz, 'vcx':vcx, 'vcy':vcy, 'vcz':vcz,'distp':distp, 'distc':distc, 'time':time}
	if filename:
		f = open(str(filename),'wb')
		pickle.dump(BHorbitInfo,f)
		f.close()
	return BHorbitInfo
			
			
	
def getBHoutput(simname,outputname,ChaNGa=True):
	getFileLists(simname)
	os.system('ls *'+outputname+'* > outfiles.list')
	files=[]
	for lines in open('outfiles.list').readlines():
        	fields=lines.split()
        	files.append(fields[0])
	print "checking for restarts..."
	if os.path.exists('restarts.txt'): os.system('rm restarts.txt')
	if len(files) > 1:
		for i in range(len(files)):
        		print i, 'file = ', files[i]
			if not ChaNGa: os.system("awk '/Restart/' "+files[i]+" >> restarts.txt")
			if ChaNGa: os.system("awk '/Restarting/' "+files[i]+" >> restarts.txt")
		ff = open('restarts.txt','r')
	if os.path.exists('out.bh'): os.system('rm out.bh')
	output = open('out.bh','a')
	for i in range(len(files)):
		print "getting BH outputs from files..."
		print i, 'file = ', files[i]
		if not ChaNGa: os.system("awk '/BHSink|Calculating/' "+files[i]+" >> tmp.bh")
		if ChaNGa: os.system("awk '/BHSink|Starting/' "+files[i]+" >> tmp.bh")
		bhf = open("tmp.bh",'r')
		if i < len(files)-1:
			restartline = ff.readline()
			if not ChaNGa: 
				N = restartline[13:-1]+'.000000'
				badline='Calculating Gravity, Step:'+N
			if ChaNGa: 
				N = restartline[14:-1]
				badline='Starting big step '+N
		else:
			badline = ''
		linarr = np.array(bhf.readlines())
		if badline: 
			print badline
			bad, = np.where(linarr == badline+'\n')
			if np.size(bad): 
				print bad, bad[0]
				linarr = linarr[np.arange(bad[0]+1)]
		bhf.close()
		os.system('rm tmp.bh')
		np.savetxt(output,linarr,fmt='%s',newline='')
	output.close()
	print "extracting data into smaller files..."
	os.system("awk '/dm/ && /dE/' out.bh > out.dm")
	os.system("awk '/C_s/' out.bh > out.cs")
	os.system("awk '/mdot/' out.bh > out.mdot")
	os.system("awk '/edible/' out.bh > out.edible")
	os.system("awk '/multi/' out.bh > out.multi")
	os.system("awk '/dist2/' out.bh > out.distance")
	os.system("awk '/CoolOffUntil/' out.bh > out.cooloff")
	os.system("awk '/Merge/' out.bh > out.merge")
	os.system("awk '/nSink/' out.bh > out.nsink")
	os.system("awk '/velocity/' out.merge > out.velocity")
	os.system("awk '/Gas|Accretion/' out.bh > out.gas")
	os.system("awk '/dx/' out.gas > out.gas.pos")
	os.system("awk '/dvx/' out.gas > out.gas.vel")
		
		
		
def read1Darray(filename,skiplines=False,dtype='int64'):
	f = open(filename,'r')
	if skiplines:
		for i in range(skiplines):
			tmp = f.readline()
	out = np.array(f.readlines()).astype(dtype)
	f.close()
	return out

#def getBHiords():
#	f = open('files.list','r')
#	files = f.readlines()
#	s = pynbody.load(files[len(files)-1].strip('\n'))
#	bhinds, = np.where(s.stars['tform']<0)
#	bhinds = len(s.dm)+len(s.gas)+bhinds
#	iord = read1Darray(files[len(files)-1].strip('\n')+'.iord',skiplines=1)
#	igasord = read1Darray(files[len(files)-1].strip('\n')+'.igasorder',skiplines=1)
#	bhiords = iord[bhinds]
#	bhigasord = igasord[bhinds]
#	f.close()
#	del(s)
#	return bhiords,bhigasord

def getBHiords(simname):
	if not os.path.exists("BHid.list"):
		print "finding IDs for all BHs that ever existed..."
		os.system("awk '{print $1}' "+simname+".orbit > BHid.list")
		f = open("BHid.list",'r')
		id = f.readlines()
		id = np.array(id)
		id = id.astype('int')
		id = np.unique(id)
		f.close()
		os.system("rm BHid.list")
		np.savetxt("BHid.list",id)
	else:
		print "previous BHid.list file found! reading it..."
		id, = readcol.readcol("BHid.list",twod=False)

	return id

def getScaleFactor(times,s):
	redshift = np.zeros(np.size(times))
	ntimes = np.size(times)
	for tt in range(np.size(times)):
		if tt%100==0:
			print tt/np.float(np.size(times)) * 100, '% done'
                def func(z):
                        return getTime(z,s) - times.in_units('Gyr')[tt]
                try: redshift[tt] = opt.newton(func,0)
		except: 
			print "ERROR did not converge", times[tt],tt
			redshift[tt] = -1
        scaleFac = 1./(1+redshift)
	return scaleFac, redshift

def getBHMergers(simname,orbitfile=None,halofile=None,outputname=None,filename=None):
	if orbitfile:
		f = open(orbitfile,'rb')
		BHorbit = pickle.load(f)
		f.close()
	if halofile:
		f2 = open(halofile,'rb')
		BHhalo = pickle.load(f2)
		f2.close()
	if orbitfile or halofile:
		if not os.path.exists(orbitfile) or not os.path.exists(halofile):
			print "ERROR: cannot fine orbit and/or halo file"
			return
	if not os.path.exists('files.list'):
                print "files.list not found.  generating list of output files..."
                getFileLists(simname)
        files = open("files.list",'r')
        f1 = files.readlines()
        s = pynbody.load(f1[0].strip('\n'))
        munits = s['mass'].units
        posunits = s['x'].units
        velunits = s['vx'].units
        potunits = s['phi'].units
        tunits = posunits/velunits
        Eunits = munits*potunits
        files.close()

	if not os.path.exists('BHmerge.txt'):
		if outputname==None:
			os.system("awk '/BHSink/ && /Merge/ && /eating/' *out* > BHmerge.txt")
		else:
			os.system("awk '/BHSink/ && /Merge/ && /eating/' *"+outputname+"* > BHmerge.txt")
	else:
		print "WARNING: BHmerge.txt already exists for this run... using that one. Please delete if you would like it to be updated."
	a,b,ID1,c,ID2,d, Time, e, f, kvel, g, h, Mratio = readcol.readcol('BHmerge.txt',twod=False)
	del(a,b,c,d,e,f,g,h)
	ID1 = np.int(ID1)
	ID2 = np.int(ID2)
	Time = np.float(Time)
	kvel = np.float(kvel)
	Mratio = np.float(Mratio)
	id2tmp, count = np.unique(ID2,return_counts=True)
	bad, = np.where(count>1)
	if len(bad)>0:
		print "Warning! Found a double counted merger. Fixing..."
		idbad = id2tmp[bad]
		for i in idbad:
			baddat, = np.where(ID2==idbad)
			np.delete(ID2,baddat[0:len(ID2)-1])
	                np.delete(ID1,baddat[0:len(ID2)-1])
	                np.delete(Time,baddat[0:len(ID2)-1])
	                np.delete(kvel,baddat[0:len(ID2)-1])
	                np.delete(Mratio,baddat[0:len(ID2)-1])
	o = np.argsort(Time)
	Time = array.SimArray(Time[o],tunits)
	Time = Time.in_units('Gyr')
	ID1 = ID1[o]
	ID2 = ID2[o]
	kvel = kvel[o]
	Mratio = Mratio[o]
	nMergers = len(Time)
	print "found", nMergers, "BH-BH mergers occuring in simulation"
	M1 = array.SimArray(np.zeros(nMergers),'Msol')
	M2 = array.SimArray(np.zeros(nMergers),'Msol')
	HaloID1 = np.zeros(nMergers)
	HaloID2 = np.zeros(nMergers)
	HaloMass1 = array.SimArray(np.zeros(nMergers),'Msol')
	HaloGas1 = array.SimArray(np.zeros(nMergers),'Msol')
	HaloStars1 = array.SimArray(np.zeros(nMergers),'Msol')
	HaloMass2 = array.SimArray(np.zeros(nMergers),'Msol')
	HaloGas2 = array.SimArray(np.zeros(nMergers),'Msol')
	HaloStars2 = array.SimArray(np.zeros(nMergers),'Msol')
	
	if BHorbit or BHhalo:
		for i in range(nMergers):
			if BHorbit:
				no1, = np.where(BHorbit['iord']==ID1[i])
				no2, = np.where(BHorbit['iord']==ID2[i])
				to, = np.where(BHorbit['data'][no1]['Time'].in_units('Gyr')<Time[i])
				if BHorbit['data'][no2]['Time'][-1].in_units('Gyr') > Time[1]:
                                	print "WARNING larger time in orbit file for BH", BHhalo['iord'][no2]," Tmerge", Time[i], "Torbit", BHorbit['data'][n1]['Time'].max()
	                        M1[i] = BHorbit['data'][no1]['mass'][to[-1]].in_units('Msol')
	                        M2[i] = BHorbit['data'][no2]['mass'][-1].in_units('Msol')
			if BHhalo:
				nh1, = np.where(BHhalo['iord']==ID1[i])
				nh2, = np.where(BHhalo['iord']==ID2[i])
				nonz, = np.where(BHhalo['mass'][nh2]>0)
				HaloID1[i] = BHhalo['haloID'][nh1][nonz[-1]]
	                        HaloID2[i] = BHhalo['haloID'][nh2][nonz[-1]]
	                        HaloMass1[i] = BHhalo['halomass'][nh1][nonz[-1]]
	                        HaloMass2[i] = BHhalo['halomass'][nh2][nonz[-1]]
	                        HaloGas1[i] = BHhalo['halogasmass'][nh1][nonz[-1]]
	                        HaloGas2i[i] = BHhalo['halogasmass'][nh2][nonz[-1]]
	                        HaloStars1[i] = BHhalo['halostarmass'][nh1][nonz[-1]]
	                        HaloStars2[i] = BHhalo['halostarmass'][nh2][nonz[-1]]
	BHmerge = {'Time':Time,'ID1':ID1,'ID2':ID2,'M1':M1,'M2':M2,'halo1':HaloID1,'halo2':HaloID2,'Hmass1':HaloMass1,'HGasMass1':HaloGas1,'HStarMass1':HaloStars1,'Hmass1':HaloMass2,'HGasMass1':HaloGas2,'HStarMass1':HaloStars2,'kickV':kvel,'ratio':Mratio}

	if filename:
		f = open(filename,'wb')
		pickle.dump(BHmerge, f)
		f.close()
	return BHmerge

def getBHorbit(simname,BHlist=[],filename=None):
	if not os.path.exists('files.list'):
		print "files.list not found.  generating list of output files..."
		getFileLists(simname)
	print "getting all BH id numbers..."
	bhids = getBHiords(simname)
	files = open("files.list",'r')
	f1 = files.readlines()
	s = pynbody.load(f1[0].strip('\n'))
	munits = s['mass'].units
	posunits = s['x'].units
	velunits = s['vx'].units
	potunits = s['phi'].units
	tunits = posunits/velunits
	Eunits = munits*potunits
	scaleUnit = pynbody.units.Unit('a')
	files.close()
	print posunits/scaleUnit
	print velunits/scaleUnit

	orbitfile = simname+".orbit"
	#print "reading "+orbitfile+"...."
	#bhorbitData = readcol.readcol(orbitfile)
	#bhids = np.unique(bhorbitData[:,0])
	if len(BHlist)>0:
                matches = np.in1d(bhids, BHlist)
                bhorbit = {'iord':bhids[matches],'data':np.array([])}
        else: bhorbit = {'iord':bhids,'data':np.array([])}
        print "there are ", len(bhids), " BHs that have existed in this simulation"
        if len(BHlist)>0: nBHs = len(BHlist)
        else: nBHs = len(bhids)

	print "getting data...."
	cnt = 0
	os.system("awk -F ' ' '{print >$1}' "+orbitfile)
	print bhorbit['iord']
	for id in bhids:
		print "getting data for BH ", id
		if len(BHlist)>0:
                        match, = np.where(BHlist==id)
                        if len(match)==0:
                                os.system("rm "+str(np.int(id)))
                                continue
		bhorbitData = readcol.readcol(str(np.int(id)))
		os.system("rm "+str(np.int(id)))
		bad, = np.where(bhorbitData[:,0] != id)
		if len(bad)>0:
			print "WARNING: bad ID found in miniorbit.txt file after awk... deleting"
			bhorbitData = np.delete(bhorbitData,bad,axis=0)
		cnt += 1
		GoodScale = True
		print "BH #"+str(cnt)+"/"+str(len(bhids))
#		curbh, = np.where(bhorbitData[:,0]==id)
		time = array.SimArray(bhorbitData[:,1],tunits)
		step = bhorbitData[:,2]
		mass = bhorbitData[:,3]
		x = bhorbitData[:,4]
		y = bhorbitData[:,5]
		z = bhorbitData[:,6]
		vx = bhorbitData[:,7]
		vy = bhorbitData[:,8]
		vz = bhorbitData[:,9]
		pot = bhorbitData[:,10]
		mdot = bhorbitData[:,11]
		deltaM = bhorbitData[:,12]
		E = bhorbitData[:,13]
		dtEff = bhorbitData[:,14]
		if len(bhorbitData[0,:])<16: 
			print "uh oh, trying to find scale factor data, but cannot!"
			scaleFac = np.ones(len(bhorbitData[:,1]))
			redshift = np.ones(len(bhorbitData[:,1]))
			GoodScale = False
		else:
			scaleFac =  bhorbitData[:,15]
			redshift = 1/scaleFac - 1
		o = np.argsort(time)
		timeOrd = time[o]
		t1 = timeOrd[0:len(timeOrd)-1]
                t2 = timeOrd[1:len(timeOrd)]
		bad = np.where(np.equal(t1,t2))
		np.delete(o,bad)
		time = array.SimArray(time[o],tunits)
		step = step[o]
		mass = array.SimArray(mass[o],munits)
		x = array.SimArray(x[o]*scaleFac[o],posunits/scaleUnit)
		y = array.SimArray(y[o]*scaleFac[o],posunits/scaleUnit)
		z = array.SimArray(z[o]*scaleFac[o],posunits/scaleUnit)
		vx = array.SimArray(vx[o]*scaleFac[o],velunits/scaleUnit)
		vy = array.SimArray(vy[o]*scaleFac[o],velunits/scaleUnit)
		vz = array.SimArray(vz[o]*scaleFac[o],velunits/scaleUnit)
		pot = array.SimArray(pot[o],potunits)
		mdot = array.SimArray(mdot[o],munits/tunits)
		deltaM = array.SimArray(deltaM[o],munits)
		E = array.SimArray(E[o],Eunits)
		dtEff = array.SimArray(dtEff[o],tunits)
		scaleFac = scaleFac[o]
		redshift = redshift[o]
		if GoodScale:
			data = {'Time':time,'step':step,'mass':mass.in_units('Msol'),'x':x.in_units('kpc'),'y':y.in_units('kpc'),'z':z.in_units('kpc'),'vx':vx.in_units('km s**-1'),'vy':vy.in_units('km s**-1'),'vz':vz.in_units('km s**-1'),'pot':pot,'mdot':mdot.in_units('Msol yr**-1'),'dM':deltaM.in_units('Msol'),'E':E,'dt':dtEff,'redshift':redshift,'scaleFac':scaleFac}
		else:
			data = {'Time':time,'step':step,'mass':mass.in_units('Msol'),'x':x,'y':y,'z':z,'vx':vx,'vy':vy,'vz':vz,'pot':pot,'mdot':mdot.in_units('Msol yr**-1'),'dM':deltaM.in_units('Msol'),'E':E,'dt':dtEff,'redshift':redshift,'scaleFac':scaleFac}
		bhorbit['data'] = np.append(bhorbit['data'],data)

	del(s)
	if filename:
                f = open(str(filename),'wb')
                pickle.dump(bhorbit,f)
                f.close()	
	return bhorbit

def LumCutOrbit(bhorbit,steplist=None,Lmax=None,Lmin=1e43,filename='bhorbitCUT.pkl'):
	if steplist:
		f = open(steplist,'r')
		steps = f.readlines()
		f.close()
		steplist = np.array(steps).astype(np.float)
	nBHs = len(bhorbit['iord'])
	bhorbitCUT = {'iord':np.array([]),'data':np.array([])}
	if not Lmax and not Lmin:
		print "both limits set to None. Not doing anything..."
		return
	cnt = 0
	for i in range(nBHs):
		if i%10 == 0: print float(i)/float(nBHs) * 100, "% done"
		if steplist: 
			match = np.in1d(bhorbit['data'][i]['step'][1:-1],steplist)
			testm, = np.where(match)
			if len(testm)==0: continue
			lum = bhorbit['data'][i]['mdot'][1:-1][testm].in_units('g s**-1')*(2.99792e10)**2*0.1
		else:
			lum = bhorbit['data'][i]['mdot'][1:-1].in_units('g s**-1')*(2.99792e10)**2*0.1
		if Lmax and Lmin: o, = np.where((lum<Lmax)&(lum>Lmin))
		if Lmax and not Lmin: o, = np.where(lum<Lmax)
		if Lmin and not Lmax: o, = np.where(lum>Lmin)
		if len(o) == 0: continue
		bhorbitCUT['iord'] = np.append(bhorbitCUT['iord'],bhorbit['iord'][i])
		bhorbitCUT['data'] = np.append(bhorbitCUT['data'],bhorbit['data'][i])
		if steplist: bhorbitCUT['data'][-1]['ind'] = testm[o]+1 
		else: bhorbitCUT['data'][-1]['ind'] = o + 1
	if filename:
		f = open(filename,'wb')
		pickle.dump(bhorbitCUT,f)
		f.close()
	return bhorbitCUT

def smoothAcc(bhorbit,trange=None,bins=100,tunits='Gyr'):
	if trange==None:
		trange = [0,bhorbit['Time'].in_units(tunits).max()]
	n,bins = np.histogram(bhorbit['Time'].in_units(tunits),bins=bins,range=trange)
	sum,bins = np.histogram(bhorbit['Time'].in_units(tunits),bins=bins,range=trange,weights=bhorbit['mdot'])
	mean = sum/n
	return mean,bins[0:len(bins)-1]+0.5*(bins[1:len(bins)]-bins[0:len(bins)-1])

def getBHFormInfo():
	f= open("files.list")
	files = f.readlines()
	s = pynbody.load(files[len(files)-1].strip())
	s.read_starlog()
	bhinds, =  np.where(s.stars['tform']<0)
	massform = s.stars['massform'][bhinds].in_units('Msol')
	tform = -1.0*s.stars['tform'][bhinds].in_units('Gyr')
	scaleFac,redshift = getScaleFactor(tform,s)
	rhoform = s.stars['rhoform'][bhinds]
	posform = s.stars['posform'][bhinds]
	tempform = s.stars['tempform'][bhinds]
	velform = s.stars['velform'][bhinds]
	forminfo = {'massform':massform, 'tform': tform,'scaleFac':scaleFac,'redform':redshift,'tempform':tempform,'velform':velform,'posform':posform,'rhoform':rhoform}
	
	return forminfo

def getBHhalo(simname,findcenter='hyb',minHM = 1e10,minNum=30,filename=None, initFile=None):
	if not os.path.exists("grpfiles.list"):
		simname_split = simname.split('.')
        	num = len(simname_split)
		os.system('ls '+simname+'.00*.grp | cut -d "." -f1-'+str(num+1)+ '> grpfiles.list' )
	if filename:
		if os.path.exists(filename):
			print "file", filename, "already exists! reading it in and appending it with new data"
			f = open(filename,'rb')
			BHhaloOLD = pickle.load(f)
			f.close()
			startStep = len(BHhaloOLD['haloID'][0])
			os.system('rm '+filename)
			print "old file has", startStep, "halos already completed"
		else:
			startStep = 0
	if initFile:
		if os.path.exists(initFile):
			print "found file ", initFile, "reading it in now"
			f = open(initFile,'rb')
                        BHhaloOLD = pickle.load(f)
                        f.close()
                        startStep = len(BHhaloOLD['haloID'][0])
                        print "given file has", startStep, "halos already completed"
			if initFile==filename:
				print "given file has same name as target file... deleting old file to replace with new file"
				os.system('rm '+filename)

	if initFile==None: startStep = 0

	f= open("grpfiles.list")
	munits = 'Msol'
	vunits = 'km s**-1'
	posunits = 'kpc'
	cposunits = 'a kpc'

	print "finding BH iords..."
        bhiords = getBHiords(simname)
	files = f.readlines()
	f.close()
	nsteps = len(files) - startStep
	nbh = len(bhiords)

	bhmass = array.SimArray(np.zeros((nbh,nsteps)),munits)
	haloid = np.zeros((nbh,nsteps))
	mhalo = array.SimArray(np.zeros((nbh,nsteps)),munits)
	mdark = array.SimArray(np.zeros((nbh,nsteps)),munits)
	mstar = array.SimArray(np.zeros((nbh,nsteps)),munits)
	mgas = array.SimArray(np.zeros((nbh,nsteps)),munits)
	#vhalo = array.SimArray(np.zeros((nbh,nsteps,3)),vunits)
	dist = array.SimArray(np.zeros((nbh,nsteps)),posunits)
	distcen = array.SimArray(np.zeros((nbh,nsteps)),posunits)
	bhpos = array.SimArray(np.zeros((nbh,nsteps,3)),posunits)
	bhposcen = array.SimArray(np.zeros((nbh,nsteps,3)),posunits)
	bhvel = array.SimArray(np.zeros((nbh,nsteps,3)),vunits)
	bhvelcen = array.SimArray(np.zeros((nbh,nsteps,3)),vunits)
	halorad = array.SimArray(np.zeros((nbh,nsteps)),posunits)
	scaleFac = np.zeros((nbh,nsteps))
	interact = np.zeros((nbh,nsteps))
	intpos = array.SimArray(np.zeros((nbh,nsteps,3)),posunits)
	intvel = array.SimArray(np.zeros((nbh,nsteps,3)),vunits)
	intdist = array.SimArray(np.zeros((nbh,nsteps)),posunits)
	#rho = array.SimArray(np.zeros((nbh,nsteps)),'g cm**-3')
	#cs = array.SimArray(np.zeros((nbh,nsteps)),'cm s**-1')


	for stepcnt in range(nsteps):
		line = files[stepcnt+startStep].strip()
		print "getting halo information for ", line
		s = pynbody.load(line)
		s.physical_units()
		cboxsize = 2*s['x'].in_units('a kpc').max()
		simBH, = np.where(np.in1d(s.star['iord'],bhiords))
		if not len(simBH):
                        print "no BHs in this step! moving on..."
                        continue
		boxsize = cboxsize.in_units('kpc')
		amigastat = readcol.readcol(line+'.amiga.stat',asdict=True)
		amigastat['cen'] = pynbody.array.SimArray((np.array([amigastat['Xc'],amigastat['Yc'],amigastat['Zc']]).T*1e3 - cboxsize/2.)*s.properties['a'],posunits)
		h = s.halos()
		#simBH, = np.where(np.in1d(s.star['iord'],bhiords))
		okgrp, = np.where(np.in1d(s.star['amiga.grp'][simBH],amigastat['Grp']))
		simBH = simBH[okgrp]
		asort = np.argsort(s.star['iord'][simBH])
		simBH = simBH[asort]
		simoutBH, = np.where(np.in1d(bhiords,s.star['iord'][simBH]))
                #outBH = np.where(np.in1d(bhiords,s.star['iord'][simBH]))
		print "there are ", len(simBH), "BHs in the step"
		allHaloID,invInd = np.unique(s.star['amiga.grp'][simBH],return_inverse=True)
		statind, = np.where(np.in1d(amigastat['Grp'],allHaloID))
		#bad, = np.where(np.in1d(allHaloID,amigastat['Grp'])==False)
		#np.delete(allHaloID,bad)
		#badind, = np.where(np.in1d(invInd,bad))
		#np.delete(invInd,badind)
		if not np.array_equal(allHaloID[invInd],amigastat['Grp'][statind[invInd]]):
			print "fuck!"
			return
		haloid[simoutBH,stepcnt] = allHaloID[invInd]
		mhalo[simoutBH,stepcnt] = pynbody.array.SimArray(amigastat['Mvir(M_sol)'][statind[invInd]],munits)
		mstar[simoutBH,stepcnt] = pynbody.array.SimArray(amigastat['StarMass(M_sol)'][statind[invInd]],munits)
		mgas[simoutBH,stepcnt]  = pynbody.array.SimArray(amigastat['GasMass(M_sol)'][statind[invInd]],munits)
		mgas[simoutBH,stepcnt]  = pynbody.array.SimArray(amigastat['GasMass(M_sol)'][statind[invInd]],munits)
		halorad[simoutBH,stepcnt] = pynbody.array.SimArray(amigastat['Rvir(kpc)'][statind[invInd]]*s.properties['a'],posunits)
		scaleFac[simoutBH,stepcnt] = s.properties['a']
		vel = np.array([amigastat['VXc'][statind[invInd]],amigastat['VYc'][statind[invInd]],amigastat['VZc'][statind[invInd]]]).T 
		bhvel[simoutBH,stepcnt,:] = s.stars['vel'][simBH].in_units(vunits) - vel 
		postemp = s.stars['pos'][simBH].in_units(posunits) - amigastat['cen'][statind[invInd]]
		postemp[(np.abs(postemp)>boxsize/2.)] = -1.0*(postemp[(np.abs(postemp)>boxsize/2.)]/np.abs(postemp[(np.abs(postemp)>boxsize/2.)])) * (boxsize-np.abs(postemp[(np.abs(postemp)>boxsize/2.)]))
		bhpos[simoutBH,stepcnt,:] = postemp
		
		bhmass[simoutBH,stepcnt] = s.stars['mass'][simBH].in_units(munits)
		dist[simoutBH,stepcnt] = np.sqrt((bhpos[simoutBH,stepcnt,:]**2).sum(axis=1))
		for cnt in range(len(allHaloID)):
			if allHaloID[cnt]==0: continue
			print allHaloID[cnt]
			oo, = np.where(amigastat['Grp']==allHaloID[cnt])
			#cen = pynbody.array.SimArray([amigastat['Xc'][oo[0]],amigastat['Yc'][oo[0]],amigastat['Zc'][oo[0]]],posunits)
			if amigastat['Mvir(M_sol)'][(amigastat['Grp']==allHaloID[cnt])]<minHM and allHaloID[cnt] > minNum: continue
			okcenter = 1
                        try:
                                pynbody.analysis.halo.center(h[allHaloID[cnt]],mode=findcenter,wrap=True,cen_size='2 kpc')
                        except ValueError:
                                okcenter = 0
                                pynbody.analysis.halo.center(h[allHaloID[cnt]],mode=findcenter,Wrap=True,cen_size='2 kpc',vel=False)	

			haloBHs, = np.where(np.in1d(h[allHaloID[cnt]].star['iord'],bhiords))
			outBH, = np.where(np.in1d(bhiords,h[allHaloID[cnt]].star['iord'][haloBHs]))
			#pynbody.transformation.inverse_translate(s,cen)
			closeBHs, = np.where((s.star['r'][simBH].in_units('kpc')<amigastat['Rvir(kpc)'][oo]*s.properties['a'])&(s.star['amiga.grp'][simBH]>allHaloID[cnt]))
			closeBHs = simBH[closeBHs]
			otheroutBHs, = np.where(np.in1d(bhiords,s.star['iord'][closeBHs]))
			
			bhposcen[outBH,stepcnt,:] = h[allHaloID[cnt]].stars['pos'][haloBHs].in_units(posunits)
			distcen[outBH,stepcnt] = h[allHaloID[cnt]].stars['r'][haloBHs].in_units(posunits)
			interact[otheroutBHs,stepcnt] = allHaloID[cnt]
			intpos[otheroutBHs,stepcnt,:] = s.stars['pos'][closeBHs].in_units(posunits)
			#intvel[otheroutBHs,stepcnt,:] = s.stars['vel'][closeBHs].in_units(vunits)
			intdist[otheroutBHs,stepcnt] = s.stars['r'][closeBHs].in_units(posunits)
			if okcenter == 1:
				intvel[otheroutBHs,stepcnt,:] = s.stars['vel'][closeBHs].in_units(vunits)
				bhvelcen[outBH,stepcnt,:] = h[allHaloID[cnt]].stars['vel'][haloBHs].in_units(vunits)

		print "deleting stuff"
		del(s)
		del(h)
		gc.collect()
		
	
	bhhalo = {'iord':bhiords,'mass':bhmass,'pos':bhpos,'poscen':bhposcen,'vel':bhvel,'velcen':bhvelcen,'haloID':haloid,'halomass': mhalo,'halostarmass':mstar,'halodarkmass':mdark,'halogasmass':mgas,'rhalo':halorad,'dist':dist,'distcen':distcen,'interact':interact,'intdist':intdist,'intvel':intvel,'intpos':intpos,'scaleFac':scaleFac}
	if startStep != 0:
		bhhalo['mass'] = np.append(BHhaloOLD['mass'], bhhalo['mass'],axis=1)
		bhhalo['pos'] = np.append(BHhaloOLD['pos'], bhhalo['pos'],axis=1)
		bhhalo['poscen'] = np.append(BHhaloOLD['poscen'], bhhalo['poscen'],axis=1)
		bhhalo['vel'] = np.append(BHhaloOLD['vel'], bhhalo['vel'],axis=1)
		bhhalo['velcen'] = np.append(BHhaloOLD['velcen'], bhhalo['velcen'],axis=1)
		bhhalo['haloID'] = np.append(BHhaloOLD['haloID'], bhhalo['haloID'],axis=1)
		bhhalo['halomass'] = np.append(BHhaloOLD['halomass'], bhhalo['halomass'],axis=1)
		bhhalo['halostarmass'] = np.append(BHhaloOLD['halostarmass'], bhhalo['halostarmass'],axis=1)
		bhhalo['halodarkmass'] = np.append(BHhaloOLD['halodarkmass'], bhhalo['halodarkmass'],axis=1)
		bhhalo['halogasmass'] = np.append(BHhaloOLD['halogasmass'], bhhalo['halogasmass'],axis=1)
		bhhalo['rhalo'] = np.append(BHhaloOLD['rhalo'], bhhalo['rhalo'],axis=1)
		bhhalo['dist'] = np.append(BHhaloOLD['dist'], bhhalo['dist'],axis=1)
		bhhalo['distcen'] = np.append(BHhaloOLD['distcen'], bhhalo['distcen'],axis=1)
		bhhalo['interact'] = np.append(BHhaloOLD['interact'],bhhalo['interact'],axis=1)
		bhhalo['intdist'] = np.append(BHhaloOLD['intdist'],bhhalo['intdist'],axis=1)
		bhhalo['intpos'] = np.append(BHhaloOLD['intpos'],bhhalo['intpos'],axis=1)
		bhhalo['intvel'] = np.append(BHhaloOLD['intvel'],bhhalo['intvel'],axis=1)
		bhhalo['scaleFac'] = np.append(BHhaloOLD['scaleFac'], bhhalo['scaleFac'],axis=1)
	if filename:
                f = open(str(filename),'wb')
                pickle.dump(bhhalo,f)
                f.close()
	return bhhalo

def getAccDens(simname,vol = 25.**3, filename='AccDens.pkl',Mlimit=1.5e6):
	munits = pynbody.units.Unit(np.str(1.9911e15)+' Msol')
        #posunits = s['x'].units
        #velunits = s['vx'].units
	#tunits = posunits/velunits
	#del(s)
	#gc.collect()
	if not os.path.exists(simname+'.BHAccLog.abridged'):#+np.str(Mlimit)):
		print "Makeing abridged Accretion log file..."
		Mlimitsim = 1.5e6/float(munits)
		cstr = """ awk '{if ($5 > """+str(Mlimitsim)+""") print $5 " " $9 " " $20}' """ + simname + ".BHAccLog > " + simname + ".BHAccLog.abridged"#+np.str(Mlimit)
		os.system(cstr)
	print "reading in data..."
	mass,dM,scale= readcol.readcol(simname+'.BHAccLog.abridged',twod=False)
	print "done!"
	#del(iord)
	#gc.collect()
	print "sorting time..."
	o = np.argsort(scale)
	#del(time)
	#gc.collect()
	print "sorting other stuff..."
	dM = pynbody.array.SimArray(dM[o],munits)
	mass = pynbody.array.SimArray(mass[o],munits)
	#time = pynbody.array.SimArray(time[o],tunits)
	scale = scale[o]
	print "summing..."
	rhoBH = np.cumsum(dM[(mass.in_units('Msol')>Mlimit)].in_units('Msol'))/vol
	scale = scale[(mass.in_units('Msol')>Mlimit)]
	del(mass)
	gc.collect()
	#time = time[(mass.in_units('Msol')>Mlimit)]
	if filename:
		print "saving data..."
		f = open(filename,'wb')
		pickle.dump([rhoBH,scale],f)
		f.close()
	return rhoBH,scale
	
def plotAccDens_v_z(rhoBH,scale,data=True,style='b-',ylog=True,xlog=True,overplot=False,lw=2,label=False):
	shankar09L = 3.2e5
	shankar09H = 5.4e5
	Salvaterra12 = 0.66e4 
	Salvaterra12zH = 9
	Salvaterra12zL = 5
	Treister13 = np.array([851.,666.,674.])
	Treister13z = np.array([6.5,7.5,8.5])
	Treister13zErr = np.array([.5,.5,.5])
	Hopkins07zp1,Hopkins07 = readcol.readcol("/u/sciteam/tremmel/QSOdata/RhoAccZ.csv",twod=False)
	Hopkins07zp1H,Hopkins07H = readcol.readcol("/u/sciteam/tremmel/QSOdata/RhoAccZPLUS.csv",twod=False)
        Hopkins07zp1L,Hopkins07L = readcol.readcol("/u/sciteam/tremmel/QSOdata/RhoAccZMINUS.csv",twod=False)
	Hopkins07perr = 10**Hopkins07H - 10**Hopkins07
	Hopkins07merr = 10**Hopkins07 - 10**Hopkins07L
	plt.plot(scale**-1,rhoBH,style,linewidth=lw,label=label)
	if data:
		shankar09 = (shankar09H +shankar09L) / 2.
		err = shankar09H - shankar09
		plt.errorbar([1.03],[shankar09],yerr=[err],color='black',fmt='D',label="Shankar+ 09")
		Salvaterra12z = (Salvaterra12zH + Salvaterra12zL)/2.
		plt.errorbar([Salvaterra12z+1],[Salvaterra12],color='black',fmt='x',xerr=[Salvaterra12zH-Salvaterra12z],yerr=0.5*Salvaterra12,lolims=[True],label='Salvaterra+ 12')
		plt.errorbar(Treister13z,Treister13,color='black',fmt='o',xerr=Treister13zErr,yerr=0.5*Treister13,lolims=[True,True,True], label='Treister+ 13')
		plt.errorbar(Hopkins07zp1,10**Hopkins07,color='grey',fmt='o',yerr=(Hopkins07merr,Hopkins07perr),label='Hopkins+ 07')
	if not overplot:
		if ylog: plt.yscale('log',base=10)
		if xlog: plt.xscale('log',base=10)
		plt.xticks([1,2,3,4,5,6,7,8,9,10],['0','1','2','3','4','5','6','7','8','9'])
	return
	
def getAccretion(simname, BHlist=[],filename=False, allData=False):
	if not os.path.exists('files.list'):
                print "files.list not found.  generating list of output files..."
                getFileLists(simname)
	bhids = getBHiords(simname)
        files = open("files.list",'r')
        f1 = files.readlines()
        s = pynbody.load(f1[0].strip('\n'))
        munits = s['mass'].units
        posunits = s['x'].units
        velunits = s['vx'].units
        potunits = s['phi'].units
        tunits = posunits/velunits
        Eunits = munits*potunits
        files.close()

	print "separating BH data..."
	acclogFile = simname+'.BHAccLog'
	os.system("awk -F ' ' '{print >$1}' "+acclogFile)
	#bhAccData = readcol.readcol(acclogFile)
	#bhids = np.unique(bhAccData[:,0])
        if len(BHlist)>0:
		matches = np.in1d(bhids, BHlist)
		bhAccHist = {'iord':bhids[matches],'data':np.array([])} 
	else: bhAccHist = {'iord':bhids,'data':np.array([])}
        print "there are ", len(bhids), " BHs that have existed in this simulation"
	if len(BHlist)>0: nBHs = len(BHlist)
	else: nBHs = len(bhids)
	print "getting data...."
        cnt = 0
        for id in bhids:
		if len(BHlist)>0:
                	match, = np.where(BHlist==id)
			if len(match)==0:
				os.system("rm "+str(np.int(id)))
				continue
		print "getting data for BH ", id
                bhAccData = readcol.readcol(str(np.int(id)))
                os.system("rm "+str(np.int(id)))
                bad, = np.where(bhAccData[:,0] != id)
                if len(bad)>0:
                        print "WARNING: bad ID found in miniorbit.txt file after awk... deleting"
                        bhAccData = np.delete(bhAccData,bad,axis=0)
                cnt += 1
                GoodScale = True
                print "BH #"+str(cnt)+"/"+str(nBHs)
		
		time = bhAccData[:,2]
		o = np.argsort(time)
                timeOrd = time[o]
                t1 = timeOrd[0:len(timeOrd)-1]
                t2 = timeOrd[1:len(timeOrd)]
                bad = np.where(np.equal(t1,t2))
                np.delete(o,bad)
                time = array.SimArray(time[o],tunits)

		iGasOrd = bhAccData[o,1]
		MgasInit = array.SimArray(bhAccData[o,3],munits)
		MbhInit = array.SimArray(bhAccData[o,4],munits)
		MgasFinal = array.SimArray(bhAccData[o,5],munits)
		MbhFinal = array.SimArray(bhAccData[o,6],munits)
		dMgas = array.SimArray(bhAccData[o,7],munits)
		dMBH = array.SimArray(bhAccData[o,8],munits)
		dMneed = array.SimArray(bhAccData[o,9],munits)
		scaleFac = bhAccData[o,19]
		dx = array.SimArray(bhAccData[o,10],posunits)*scaleFac
		dy = array.SimArray(bhAccData[o,11],posunits)*scaleFac
		dz = array.SimArray(bhAccData[o,12],posunits)*scaleFac
		dvx = array.SimArray(bhAccData[o,13],velunits)
                dvy = array.SimArray(bhAccData[o,14],velunits)
                dvz = array.SimArray(bhAccData[o,15],velunits)
		Ugas = array.SimArray(bhAccData[o,16],Eunits)
		fBall = array.SimArray(bhAccData[o,17],posunits)*scaleFac
		tCoolOff = array.SimArray(bhAccData[o,18],tunits)
		density = array.SimArray(bhAccData[o,20],munits/posunits**3)*scaleFac**(-3)
		temp = array.SimArray(bhAccData[o,21],'K')
		metals = array.SimArray(bhAccData[o,22])
		if allData:
			datastruct = {'time':time.in_units('Gyr'),'Mgas':MgasInit.in_units('Msol'),'Mbh':MbhInit.in_units('Msol'),'MgasFinal':MgasFinal.in_units('Msol'),'MbhFinal':MbhFinal.in_units('Msol'),'deltaMgas':dMgas.in_units('Msol'),'deltaM':dMBH.in_units('Msol'),'Mneed':dMneed.in_units('Msol'),'dx':dx.in_units('kpc'),'dy':dy.in_units('kpc'),'dz':dz.in_units('kpc'),'dvx':dvx.in_units('kpc'),'dvy':dvy.in_units('kpc'),'dvz':dvz.in_units('kpc'),'Ugas':Ugas,'fBall':fBall.in_units('kpc',a=1),'tCoolOff':tCoolOff,'scaleFac':scaleFac,'density':density.in_units('m_p cm**-3',a=1),'temp':temp,'metals':metals}
		else:
			datastruct =  {'time':time.in_units('Gyr'),'Mgas':MgasInit.in_units('Msol'),'Mbh':MbhInit.in_units('Msol'),'deltaM':dMBH.in_units('Msol'),'dx':dx.in_units('kpc',a=1),'dy':dy.in_units('kpc',a=1),'dz':dz.in_units('kpc',a=1),'dvx':dvx.in_units('km s**-1',a=1),'dvy':dvy.in_units('km s**-1',a=1),'dvz':dvz.in_units('km s**-1',a=1),'Ugas':Ugas,'fBall':fBall.in_units('kpc',a=1),'tCoolOff':tCoolOff,'scaleFac':scaleFac,'density':density.in_units('m_p cm**-3',a=1),'temp':temp,'metals':metals}
		bhAccHist['data'] = np.append(bhAccHist['data'],datastruct)
	del(s)		
        if filename:
                f = open(str(filename),'wb')
                pickle.dump(bhAccHist,f)
                f.close()
        return bhAccHist
	
#def BHMassDen(bhAccHist,simname):
#	sl = pynbody.tipsy.StarLog(simname+'.starlog')
#	tformBH = -1.0*sl['tform'][(sl['tform']<0)].in_units('yr')

def plotOccFrac(sim,centrals=True,rlim=1,cum=True,bins=10):
	'''
	plot the black hole occupation fraction of halos as a function of mass

	----inputs-----
	sim = name of the snapshot you wish to analyze
	centrals = whether or not you only want to count black holes within rlim from center of halo
	rlim = the maximum radius from the halo center that you would define a black hole to be "central"
	cum = whether you want a cumulative distribution
	bins = number of log bins in halo mass you want
	----outputs----
	array[occupation fraction]
	array[mass bins]
	also a nice plot!
	'''
	stat = readcol.readcol(sim+'.amiga.stat',skipline=1)
	Mvir = stat[:,5].astype('float')
	s = pynbody.load(sim)
	h = s.halos()
	bhgrps = s.stars['amiga.grp'][(s.stars['tform']<0)]
	print "calculating Min and Max halo masses..."
	MvirMin = h[bhgrps.max()]['mass'].in_units('Msol').sum()
	MvirMax = h[1]['mass'].in_units('Msol').sum()
	print "Min: ", MvirMin, " Max: ", MvirMax
	dLogM = (np.log10(MvirMax) - np.log10(MvirMin))/bins
	BHFlag = np.zeros(bhgrps.max())
	ugrps = np.unique(bhgrps)
	print "there are ", len(ugrps), "halos with BHs"
	print "determining existence of central BHs..."
	for i in ugrps:
		if i == 0: continue
		print "halo ", i
		cen = halo.shrink_sphere_center(h[i])
		h[i].stars['pos'] -= cen
		if len(h[i].stars[((h[i].stars['tform']<0)&(h[i].stars['r'].in_units('kpc')<rlim))])>0: BHFlag[i-1] = 1
		h[i].stars['pos'] += cen

	occNum, Mbins = np.histogram(np.log10(Mvir[np.arange(bhgrps.max())]),bins=bins,weights = BHFlag)
	HNum, Mbins2 = np.histogram(np.log10(Mvir[np.arange(bhgrps.max())]),bins=bins)
	if cum==True:
		occFrac = np.cumsum(occNum)/np.cumsum(HNum).astype('float')
	else:
		occFrac = occNum/HNum.astype('float')
	return Mbins,occFrac
	


class BH(object):
        def __init__(self,simname,outputname,findcenter=True,filename=False,ChaNGa=True):
                if not filename: filename = simname+'.BHinfo.pkl'
		if os.path.exists(filename):
			print "hooray, pickle files exists!"
                        ff = open(simname+'.BHinfo.pkl','r')
                        tmp = pickle.load(ff)
                        self.__dict__.update(tmp.__dict__)
			ff.close()
                else:
			print "pickle file not found... making new object"
			if not os.path.exists("out.dm"):
				print "output files not found... reading in output now..."
				getBHoutput(simname,outputname,ChaNGa=ChaNGa)
			print "pickle file not found... making new object"
			print "loading in iords, igasords for bh particles..."
			#bhiords,bhigasords = getBHiords()
			#self.iords = {'iords':bhiords,'igasords':bhigasords}
			bhigasords = 0
			print "getting orbit data..."
			self.orbit =  getBHorbit(simname,outputname)
			print "getting accretion data..."
			self.acc = getAccretion(simname)
			print "getting halo data..."
			self.halo = getBHhalo(simname,self.orbit['iord'],findcenter=findcenter) 
			print "getting formation data..."
			self.form = getBHFormInfo()
			print "saving BH structure to: ", filename
			ff = open(filename,'wb')
			pickle.dump(self, ff)
			ff.close()
		
#	x = raw_input("question:")
