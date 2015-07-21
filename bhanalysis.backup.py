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
	os.system('ls  *.den | cut -d"." -f1-'+str(num+1)+' > files.list')
        os.system('ls *.den | cut -d"." -f'+str(num+1)+' > steps.list')
	
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

def getBHiords():
	f = open('files.list','r')
	bhiords = np.array([])
	for line in f:
		print "finding iords in file ",line 
		s = pynbody.load(line.strip('\n'))
		if len(s.star)==0: continue
		bhiords = np.append(bhiords,s.stars['iord'][(s.stars['tform']<0)])
		del(s)
	bhiords = np.unique(bhiords)
	return bhiords

def getScaleFactor(times,s):
	redshift = np.zeros(np.size(times))
	for tt in range(np.size(times)):
                def func(z):
                        return getTime(z,s) - times.in_units('Gyr')[tt]
                redshift[tt] = opt.newton(func,0)
        scaleFac = 1./(1+redshift)
	return scaleFac, redshift

def getBHMergers(simname,orbitfile,halofile,outputname=None,filename=None):
	f = open(orbitfile,'rb')
	BHorbit = pickle.load(f)
	f.close()
	f2 = open(halofile,'rb')
	BHhalo = pickle.load(f2)
	f2.close()
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
		print "BHmerge.txt already exists for this run... please delete before continuing"
		return
	a,b,ID1,c,ID2,d, Time, e, f, kvel, g, h, Mratio = readcol.readcol('BHmerge.txt',twod=False)
	del(ai,b,c,d,e,f,g,h)
	o = np.artsort(Time)
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
	haloID1 = np.zeros(nMergers)
	haloID2 = np.zeros(nMergers)
	HaloMass1 = array.SimArray(np.zeros(nMergers),'Msol')
	HaloGas1 = array.SimArray(np.zeros(nMergers),'Msol')
	HaloStars1 = array.SimArray(np.zeros(nMergers),'Msol')
	HaloMass2 = array.SimArray(np.zeros(nMergers),'Msol')
        HaloGas2 = array.SimArray(np.zeros(nMergers),'Msol')
        HaloStars2 = array.SimArray(np.zeros(nMergers),'Msol')
	
	for i in range(nMergers):
		no1, = np.where(BHorbit['iord']==ID1[i])
		nh1, = np.where(BHhalo['iord']==ID1[i])
		no2, = np.where(BHorbit['iord']==ID2[i])
                nh2, = np.where(BHhalo['iord']==ID2[i])
		


		to, = np.where(BHorbit['data'][no1]['Time'].in_units('Gyr')<Time[i])
		
		if BHorbit['data'][no2]['Time'][-1].in_units('Gyr') > Time[1]: 
			print "WARNING larger time in orbit file for BH", BHhalo['iord'][no2]," Tmerge", Time[i], "Torbit", BHorbit['data'][n1]['Time'].max()
		M1[i] = BHorbit['data'][no1]['mass'][to[-1]].in_units('Msol')
		M2[i] = BHorbit['data'][no2]['mass'][-1].in_units('Msol')

		o, = np.where(BHhalo['mass'][nh2]>0)
		HaloID1[i] = BHhalo['haloID'][nh1][o[-1]]
		HaloID2[i] = BHhalo['haloID'][nh2][o[-1]]
		HaloMass1[i] = BHhalo['halomass'][nh1][o[-1]]
		HaloMass2[i] = BHhalo['halomass'][nh2][o[-1]]
		HaloGas1[i] = BHhalo['halogasmass'][nh1][o[-1]]
                HaloGas2i[i] = BHhalo['halogasmass'][nh2][o[-1]]
		HaloStar1[i] = BHhalo['halostarmass'][nh1][o[-1]]
                HaloStar2[i] = BHhalo['halostarmass'][nh2][o[-1]]
	BHmerge = {'Time':Time,'M1':M1,'M2':M2,'halo1':HaloID1,'halo2':HaloID2,'Hmass1':HaloMass1,'HGasMass1':HaloGas1,'HStarMass1':HaloStar1,'Hmass1':HaloMass2,'HGasMass1':HaloGas2,'HStarMass1':HaloStar2,'kickV':kvel,'ratio':Mratio}

	if filename:
		f = open(filename,'wb')
		pickle.dump(BHmerge, f)
		f.close()
	return BHmerge

def getBHorbit(simname,filename=None):
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

	orbitfile = simname+".orbit"
	print "reading "+orbitfile+"...."
	bhorbitData = readcol.readcol(orbitfile)
	print "collecting IDs..."
	bhids = np.unique(bhorbitData[:,0])
	bhorbit = {'iord':bhids,'data':np.array([])}
	print "there are ", len(bhids), " BHs that have existed in this simulation"
	print "getting data...."
	cnt = 0
	for id in bhids:
		cnt += 1
		GoodScale = True
		print "BH #"+str(cnt)+"/"+str(len(bhids))
		curbh, = np.where(bhorbitData[:,0]==id)
		time = array.SimArray(bhorbitData[curbh,1],tunits)
		step = bhorbitData[curbh,2]
		mass = bhorbitData[curbh,3]
		x = bhorbitData[curbh,4]
		y = bhorbitData[curbh,5]
		z = bhorbitData[curbh,6]
		vx = bhorbitData[curbh,7]
		vy = bhorbitData[curbh,8]
		vz = bhorbitData[curbh,9]
		pot = bhorbitData[curbh,10]
		mdot = bhorbitData[curbh,11]
		deltaM = bhorbitData[curbh,12]
		E = bhorbitData[curbh,13]
		dtEff = bhorbitData[curbh,14]
		if len(bhorbitData[0,:])<16: 
			print "uh oh, trying to find scale factor data, but cannot!"
			scaleFac = np.zeros(len(curbh))
			redshift = np.zeros(len(curbh))
			GoodScale = False
		else:
			scaleFac =  bhorbitData[curbh,15]
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
		x = array.SimArray(x[o],posunits)
		y = array.SimArray(y[o],posunits)
		z = array.SimArray(z[o],posunits)
		vx = array.SimArray(vx[o],velunits)
		vy = array.SimArray(vy[o],velunits)
		vz = array.SimArray(vz[o],velunits)
		pot = array.SimArray(pot[o],potunits)
		mdot = array.SimArray(mdot[o],munits/tunits)
		deltaM = array.SimArray(deltaM[o],munits)
		E = array.SimArray(E[o],Eunits)
		dtEff = array.SimArray(dtEff[o],tunits)
		scaleFac = scaleFac[o]
		redshift = redshift[o]
		if GoodScale:
			data = {'Time':time,'step':step,'mass':mass.in_units('Msol yr**-1'),'x':x.in_units('kpc',a=scaleFac),'y':y.in_units('kpc',a=scaleFac),'z':z.in_units('kpc',a=scaleFac),'vx':vx.in_units('km s**-1',a=scaleFac),'vy':vy.in_units('km s**-1',a=scaleFac),'vz':vz.in_units('km s**-1',a=scaleFac),'pot':pot,'mdot':mdot.in_units('Msol yr**-1'),'dM':deltaM.in_units('Msol'),'E':E.in_units('ergs'),'dt':dtEff,'redshift':redshift,'scaleFac':scaleFac}
		else:
			data = {'Time':time,'step':step,'mass':mass.in_units('Msol'),'x':x,'y':y,'z':z,'vx':vx,'vy':vy,'vz':vz,'pot':pot,'mdot':mdot.in_units('Msol yr**-1'),'dM':deltaM.in_units('Msol'),'E':E,'dt':dtEff,'redshift':redshift,'scaleFac':scaleFac}
		bhorbit['data'] = np.append(bhorbit['data'],data)

	del(s)
	if filename:
                f = open(str(filename),'wb')
                pickle.dump(bhorbit,f)
                f.close()	
	return bhorbit

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

def getBHhalo(simname,findcenter='mass',filename=None, initFile=None):
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

	if findcenter != 'mass' and findcenter != 'pot':
		print 'Warning: findcenter = mass or pot. Cannot understand given input, reverting to shrink_sphere_center for halo center calculation.'
		findcenter = 'mass'

	print "finding BH iords..."
	bhorbitdata = readcol.readcol(simname+'.orbit')
        bhiords = np.unique(bhorbitdata[:,0])
	files = f.readlines()
	nsteps = len(files) - startStep
	nbh = len(bhiords)
	print "cleaning up... and initializing arrays"
	del(bhorbitdata)

	bhmass = array.SimArray(np.zeros((nbh,nsteps)),munits)
	haloid = np.zeros((nbh,nsteps))
	mhalo = array.SimArray(np.zeros((nbh,nsteps)),munits)
	mdark = array.SimArray(np.zeros((nbh,nsteps)),munits)
	mstar = array.SimArray(np.zeros((nbh,nsteps)),munits)
	mgas = array.SimArray(np.zeros((nbh,nsteps)),munits)
	vhalo = array.SimArray(np.zeros((nbh,nsteps,3)),vunits)
	dist = array.SimArray(np.zeros((nbh,nsteps)),posunits)
	bhpos = array.SimArray(np.zeros((nbh,nsteps,3)),posunits)
	bhvel = array.SimArray(np.zeros((nbh,nsteps,3)),vunits)
	halocen = array.SimArray(np.zeros((nbh,nsteps,3)),posunits)
	halorad = array.SimArray(np.zeros((nbh,nsteps)),posunits)
	scaleFac = np.zeros((nbh,nsteps))
	rho = array.SimArray(np.zeros((nbh,nsteps)),'g cm**-3')
	cs = array.SimArray(np.zeros((nbh,nsteps)),'cm s**-1')


	for stepcnt in range(nsteps):
		line = files[stepcnt+startStep].strip()
		print "getting halo information for ", line
		s = pynbody.load(line)
		h = s.halos()
		if not np.size(s.star[(s.star['tform']<0)]):
			print "no BHs in this step! moving on..."
                        continue
		else: print "there are ", np.size(s.star[(s.star['tform']<0)]), "BHs in the step"
		print "finding halo centers..."
		allHaloID = np.unique(s.star['amiga.grp'][(s.star['tform']<0)])
		centers = np.zeros((len(allHaloID),3))
		velocities = np.zeros((len(allHaloID),3))
		radius = np.zeros(len(allHaloID))
		for cnt in range(len(allHaloID)):
			if findcenter=='mass': cen = halo.shrink_sphere_center(h[allHaloID[cnt]])
			if findcenter=='phi': cen = h[allHaloID[cnt]]['pos'][(h[allHaloID[cnt]]['phi']==float(h[allHaloID[cnt]]['phi'].min()))]
			vel = halo.center_of_mass_velocity(h[allHaloID[cnt]])
			h[allHaloID[cnt]]['pos'] -= cen
			rad = h[allHaloID[cnt]]['r'].max()
			centers[cnt,:] = cen.in_units(posunits,a=s.properties['a'])
			velocities[cnt,:] = vel.in_units(vunits,a=s.properties['a'])
			radius[cnt] = rad.in_units(posunits,a=s.properties['a'])
			h[allHaloID[cnt]]['pos'] += cen	

		print "getting BH data for this step..."
		for cnt in range(nbh):
			curbh, = np.where(s.star['iord']==bhiords[cnt])
			if len(curbh)==0: continue
			bhmass[cnt,stepcnt] = s.star['mass'].in_units('Msol')[curbh[0]]

			haloid[cnt,stepcnt] = s.star['amiga.grp'][curbh]
			mhalo[cnt,stepcnt] = h[haloid[cnt,stepcnt]]['mass'].sum().in_units(munits)
			mstar[cnt,stepcnt] = h[haloid[cnt,stepcnt]].star['mass'].sum().in_units(munits)
			mgas[cnt,stepcnt] = h[haloid[cnt,stepcnt]].gas['mass'].sum().in_units(munits)

			vhalo[cnt,stepcnt,:] = velocities[(allHaloID==haloid[cnt,stepcnt]),:]
			bhvel[cnt,stepcnt,:] = s.star['vel'].in_units(vunits)[curbh[0]] - vhalo[cnt,stepcnt,:]
			halocen[cnt,stepcnt,:] = centers[(allHaloID==haloid[cnt,stepcnt]),:]
			halorad[cnt,stepcnt] = radius[(allHaloID==haloid[cnt,stepcnt])]
			bhpos[cnt,stepcnt,:] = s.star['pos'].in_units(posunits)[curbh[0]]-halocen[cnt,stepcnt,:] 
			dist[cnt,stepcnt] = np.sqrt((bhpos[cnt,stepcnt,:]**2).sum())

			if len(h[haloid[cnt,stepcnt]].gas) > 0:
				h[haloid[cnt,stepcnt]].gas['pos'] -= s.star['pos'][curbh[0]]
				o = np.argsort(h[haloid[cnt,stepcnt]].gas['r'])
				rho[cnt,stepcnt] = h[haloid[cnt,stepcnt]].gas['rho'][o[0:32]].in_units('g cm**-3').mean()	
				cs[cnt,stepcnt] = h[haloid[cnt,stepcnt]].gas['cs'][o[0:32]].in_units('cm s**-1').mean()		
				h[haloid[cnt,stepcnt]].gas['pos'] += s.star['pos'][curbh[0]]

			else:
				rho[cnt,stepcnt] = 0
				cs[cnt,stepcnt] = 0

			scaleFac[cnt,stepcnt] = s.properties['a']	

		print "deleting stuff"
		del(s)
		del(h)
		
	if startStep == 0:
		bhhalo = {'iord':bhiords,'mass':bhmass,'pos':bhpos,'vel':bhvel,'haloID':haloid,'halomass': mhalo,'haloradius':halorad,'halostarmass':mstar,'halodarkmass':mdark,'halogasmass':mgas,'halocen':halocen,'halovel':vhalo,'dist':dist,'scaleFac':scaleFac,'rho':rho,'cs':cs}
	else:
		bhhalo = {'iord':bhiords,'mass':bhmass,'pos':bhpos,'vel':bhvel,'haloID':haloid,'halomass': mhalo,'haloradius':halorad,'halostarmass':mstar,'halodarkmass':mdark,'halogasmass':mgas,'halocen':halocen,'halovel':vhalo,'dist':dist,'scaleFac':scaleFac,'rho':rho,'cs':cs}
		bhhalo['mass'] = np.append(BHhaloOLD['mass'], bhhalo['mass'],axis=1)
		bhhalo['pos'] = np.append(BHhaloOLD['pos'], bhhalo['pos'],axis=1)
		bhhalo['vel'] = np.append(BHhaloOLD['vel'], bhhalo['vel'],axis=1)
		bhhalo['haloID'] = np.append(BHhaloOLD['haloID'], bhhalo['haloID'],axis=1)
		bhhalo['halomass'] = np.append(BHhaloOLD['halomass'], bhhalo['halomass'],axis=1)
		bhhalo['haloradius'] = np.append(BHhaloOLD['haloradius'], bhhalo['haloradius'],axis=1)
		bhhalo['halostarmass'] = np.append(BHhaloOLD['halostarmass'], bhhalo['halostarmass'],axis=1)
		bhhalo['halodarkmass'] = np.append(BHhaloOLD['halodarkmass'], bhhalo['halodarkmass'],axis=1)
		bhhalo['halogasmass'] = np.append(BHhaloOLD['halogasmass'], bhhalo['halogasmass'],axis=1)
		bhhalo['halocen'] = np.append(BHhaloOLD['halocen'], bhhalo['halocen'],axis=1)
		bhhalo['halovel'] = np.append(BHhaloOLD['halovel'], bhhalo['halovel'],axis=1)
		bhhalo['dist'] = np.append(BHhaloOLD['dist'], bhhalo['dist'],axis=1)
		bhhalo['scaleFac'] = np.append(BHhaloOLD['scaleFac'], bhhalo['scaleFac'],axis=1)
		bhhalo['rho'] = np.append(BHhaloOLD['rho'], bhhalo['rho'],axis=1)
		bhhalo['cs'] = np.append(BHhaloOLD['cs'], bhhalo['cs'],axis=1)
	if filename:
                f = open(str(filename),'wb')
                pickle.dump(bhhalo,f)
                f.close()
	return bhhalo

def getAccretion(simname, filename=False, getGasInfo=False, allData=False):
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

	print "reading in data..."
	acclogFile = simname+'.BHAccLog'
	bhAccData = readcol.readcol(acclogFile)
	bhids = np.unique(bhAccData[:,0])
        bhAccHist = {'iord':bhids,'data':np.array([])}
        print "there are ", len(bhids), " BHs that have existed in this simulation"
	print "getting data...."
        cnt = 0
        for id in bhids:
		cnt += 1
		print "BH #"+str(cnt)+"/"+str(len(bhids))
		curbh, = np.where(bhAccData[:,0]==id)
		
		time = bhAccData[curbh,2]
		o = np.argsort(time)
                timeOrd = time[o]
                t1 = timeOrd[0:len(timeOrd)-1]
                t2 = timeOrd[1:len(timeOrd)]
                bad = np.where(np.equal(t1,t2))
                np.delete(o,bad)
                time = array.SimArray(time[o],tunits)

		iGasOrd = bhAccData[curbh[o],1]
		MgasInit = array.SimArray(bhAccData[curbh[o],3],munits)
		MbhInit = array.SimArray(bhAccData[curbh[o],4],munits)
		MgasFinal = array.SimArray(bhAccData[curbh[o],5],munits)
		MbhFinal = array.SimArray(bhAccData[curbh[o],6],munits)
		dMgas = array.SimArray(bhAccData[curbh[o],7],munits)
		dMBH = array.SimArray(bhAccData[curbh[o],8],munits)
		dMneed = array.SimArray(bhAccData[curbh[o],9],munits)
		dx = array.SimArray(bhAccData[curbh[o],10],posunits)
		dy = array.SimArray(bhAccData[curbh[o],11],posunits)
		dz = array.SimArray(bhAccData[curbh[o],12],posunits)
		dvx = array.SimArray(bhAccData[curbh[o],13],velunits)
                dvy = array.SimArray(bhAccData[curbh[o],14],velunits)
                dvz = array.SimArray(bhAccData[curbh[o],15],velunits)
		Ugas = array.SimArray(bhAccData[curbh[o],16],Eunits)
		fBall = array.SimArray(bhAccData[curbh[o],17],posunits)
		tCoolOff = array.SimArray(bhAccData[curbh[o],18],tunits)
		scaleFac = bhAccData[curbh[o],19]
		
		if allData:
			datastruct = {'time':time,'Mgas':MgasInit,'Mbh':MbhInit,'MgasFinal':MgasFinal,'MbhFinal':MbhFinal,'deltaMgas':dMgas,'deltaM':dMBH,'Mneed':dMneed,'dx':dx,'dy':dy,'dz':dz,'dvx':dvx,'dvy':dvy,'dvz':dvz,'Ugas':Ugas,'fBall':fBall,'tCoolOff':tCoolOff,'scaleFac':scaleFac}
		else:
			datastruct =  {'time':time,'Mgas':MgasInit,'Mbh':MbhInit,'deltaM':dMBH,'dx':dx,'dy':dy,'dz':dz,'dvx':dvx,'dvy':dvy,'dvz':dvz,'Ugas':Ugas,'fBall':fBall,'tCoolOff':tCoolOff,'scaleFac':scaleFac}
		bhAccHist['data'] = np.append(bhAccHist['data'],datastruct)
	del(s)		
        if filename:
                f = open(str(filename),'wb')
                pickle.dump(bhAccHist,f)
                f.close()
        return bhAccHist

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
