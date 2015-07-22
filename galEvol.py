import volumeAnalysis
import bhanalysis
import matplotlib.colors as pltcolors
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pynbody
import matplotlib.patches as mpatches
import gc

plt.ion()
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
plt.rc('font', weight='medium')
plt.rc('axes', linewidth=2)
plt.rc('xtick.major',width=2)
plt.rc('ytick.major',width=2)

def colorcolor(halos=[1,2],average=False,simdirs=['romulus8.256gst3.bwBH'],dodust=True,dustfile = 'dust.pkl',magfile = 'mags.pkl',hcol=['red','blue'],hcmap=['Reds','Blues'],simmark=['o'],msize=100,plotData=True,simlabels=['Romulus'],hlabels=['halo 1','halo 2'],overplot=False,cbar=True):
	simcnt = 0
	for dir in simdirs:
		print "getting data for ", dir
		magf = open(dir+'/'+magfile,'rb')
		mags = pickle.load(magf)
		magf.close()
		zmax = mags['z'].max()
		zmin = mags['z'].min()
		red = mags['z']
		print zmax, zmin
		redcNorm = pltcolors.Normalize(1./zmax,1./zmin)
		if dodust:
			dustf = open(dir+'/'+dustfile,'rb')
			dust = pickle.load(dustf)
			dustf.close()
		hcnt = 0
		color1 = np.zeros((len(halos),len(mags['v'][0])))
		color2 = np.zeros((len(halos),len(mags['v'][0])))
		for h in halos:
			u, = np.where(np.array(mags['halos'])==h)
			if len(u) > 1: 
				print "WARNING more than one entry in magnitudes file for halo", h, "using first entry..."
			if len(u) == 0:
				print "Halo", h, " not found in magnitudes file... skipping..."
				continue
			color1[hcnt,:] = mags['v'][u[0]] - mags['j'][u[0]]
			color2[hcnt,:] = mags['u'][u[0]] - mags['v'][u[0]]
			if dodust==True:
				ud, = np.where(np.array(dust['halos'])==h)
				if len(u) > 1: print "WARNING more than one entry in dust file for halo", h, "using first entry..."
        	                if len(u) == 0:
	                                print "Halo", h, " not found in dust file... skipping..."
                	                continue
				reddening1 = dust['v'][ud[0]] - dust['j'][ud[0]]
				reddening2 =  dust['u'][ud[0]] - dust['v'][ud[0]]
				color1[hcnt,:] += reddening1
				color2[hcnt,:] += reddening2
#			print len(color1), len(color2), len(red), len(1./red)
#			print hcol[hcnt],hcmap[hcnt],1/red
			if average==False:
				plt.scatter(color1[hcnt,:],color2[hcnt,:],c=1./red,norm=redcNorm,cmap=hcmap[hcnt],marker=simmark[simcnt],s=msize,label=simlabels[simcnt]+" "+hlabels[hcnt],color=hcol[hcnt])
			hcnt += 1
		if average==True:
			plt.scatter(color1.mean(axis=0),color2.mean(axis=0),c=1./red,norm=redcNorm,cmap=hcmap[0],marker=simmark[simcnt],s=msize,label=simlabels[simcnt]+" "+hlabels[0],color=hcol[0])
			plt.errorbar(color1.mean(axis=0),color2.mean(axis=0),xerr=color1.std(axis=0),yerr=color2.std(axis=0),fmt='o',color=hcol[0],markersize=0,linewidth=2,elinewidth=0.75)
		simcnt += 1


	if plotData==True:
		plt.errorbar(volumeAnalysis.CANDELS_MW['VJ'],volumeAnalysis.CANDELS_MW['UV'],xerr=[volumeAnalysis.CANDELS_MW['VJ-'],volumeAnalysis.CANDELS_MW['VJ+']],yerr=[volumeAnalysis.CANDELS_MW['UV-'],volumeAnalysis.CANDELS_MW['UV+']],fmt='o-',color='grey',markersize=0,elinewidth=.75,linewidth=2)
	        plt.errorbar(volumeAnalysis.CANDELS_M31['VJ'],volumeAnalysis.CANDELS_M31['UV'],xerr=[volumeAnalysis.CANDELS_M31['VJ-'],volumeAnalysis.CANDELS_M31['VJ+']],yerr=[volumeAnalysis.CANDELS_M31['UV-'],volumeAnalysis.CANDELS_M31['UV+']],fmt='o-',color='k',markersize=0,elinewidth=.75,linewidth=2)
	        plt.scatter(volumeAnalysis.CANDELS_MW['VJ'],volumeAnalysis.CANDELS_MW['UV'],c=1./np.array(volumeAnalysis.CANDELS_MW['redshift']),norm=redcNorm,cmap='Greys',s=150,marker='D',label='CANDELS MW',linewidth=1.5,color='k')
	        plt.scatter(volumeAnalysis.CANDELS_M31['VJ'],volumeAnalysis.CANDELS_M31['UV'],c=1./np.array(volumeAnalysis.CANDELS_M31['redshift']),norm=redcNorm,cmap='Greys',s=150,marker='^',label='CANDELS M31',linewidth=1.5,color='k')
	if cbar==True:
        	cbar = plt.colorbar(ticks=[0.25,0.5,1,1.5,2.0])
        	cbar.set_label('Redshift',fontsize=30)
		cbar.set_ticklabels(['4','2','1','0.667','0.5'])

	if overplot==False:
		plt.plot([1.5,1.5],[1.9,2.7],'k--',linewidth=1.5)
	        plt.plot([1.0,1.5],[1.3,1.9],'k--',linewidth=1.5)
	        plt.plot([-0.5,1.0],[1.3,1.3],'k--',linewidth=1.5)
	        plt.xlim(-0.5,2.1)
	        plt.ylim(-0.2,2.3)
	        plt.ylabel('U-V color',fontsize=40)
	        plt.xlabel('V-J color',fontsize=40)
	        plt.xticks(fontsize=20)
	plt.legend(loc='upper left',fontsize=30)
	return

def colortime(halos=[1,2],average=False,simdirs=['romulus8.256gst3.bwBH'],dodust=True,dustfile = 'dust.pkl',magfile = 'mags.pkl',hcol=['red','blue'],simstyle=['solid','solid'],plotData=True,simlabels=['Romulus'],hlabels=['halo 1','halo 2'],overplot=False):
	simcnt = 0
        for dir in simdirs:
                print "getting data for ", dir
                magf = open(dir+'/'+magfile,'rb')
                mags = pickle.load(magf)
                magf.close()
                zmax = mags['z'].max()
                zmin = mags['z'].min()
                red = mags['z']
		if dodust:
                        dustf = open(dir+'/'+dustfile,'rb')
                        dust = pickle.load(dustf)
                        dustf.close()
                hcnt = 0
		f = open(dir+'/files.list')
		files = f.readlines()
		f.close()
		s = pynbody.load(dir+'/'+files[0].strip('\n'))
		time = [bhanalysis.getTime(z,s) for z in red]
		color = np.zeros((len(halos),len(mags['v'][0])))
                for h in halos:
                        u, = np.where(np.array(mags['halos'])==h)
                        if len(u) > 1:
                                print "WARNING more than one entry in magnitudes file for halo", h, "using first entry..."
                        if len(u) == 0:
                                print "Halo", h, " not found in magnitudes file... skipping..."
                                continue
			if average==True:
				color[hcnt,:] = mags['u'][u[0]] - mags['v'][u[0]]
			else:
				color[hcnt,:] = mags['u'][u[0]] - mags['v'][u[0]]
			if dodust==True:
                                ud, = np.where(np.array(dust['halos'])==h)
                                if len(u) > 1: print "WARNING more than one entry in dust file for halo", h, "using first entry..."
                                if len(u) == 0:
                                        print "Halo", h, " not found in dust file... skipping..."
                                        continue
                                reddening =  dust['u'][ud[0]] - dust['v'][ud[0]]
                                color += reddening
			if average==False:
				plt.plot(time,color[hcnt,:],color=hcol[hcnt],linewidth=2,linestyle=simstyle[simcnt],label=simlabels[simcnt]+" "+hlabels[hcnt])
			hcnt += 1
		if average==True:
			plt.errorbar(time,color.mean(axis=0),color=hcol[0],yerr=color.std(axis=0),linewidth=2,linestyle=simstyle[simcnt],label=simlabels[simcnt]+" "+hlabels[0])
		simcnt += 1
	
	if plotData == True:
		volumeAnalysis.CANDELS_M31['Time'] = [bhanalysis.getTime(z,s) for z in volumeAnalysis.CANDELS_M31['redshift']]
	        volumeAnalysis.CANDELS_MW['Time'] = [bhanalysis.getTime(z,s) for z in volumeAnalysis.CANDELS_MW['redshift']]
		plt.errorbar(volumeAnalysis.CANDELS_M31['Time'],volumeAnalysis.CANDELS_M31['UV'],yerr=[volumeAnalysis.CANDELS_M31['UV-'],volumeAnalysis.CANDELS_M31['UV+']],fmt='^',color='k',markersize=10,label='CANDELS M31')
	        plt.errorbar(volumeAnalysis.CANDELS_MW['Time'],volumeAnalysis.CANDELS_MW['UV'],yerr=[volumeAnalysis.CANDELS_MW['UV-'],volumeAnalysis.CANDELS_MW['UV+']],fmt='D',color='grey',markersize=10,label='CANDELS MW')
	if overplot==False:
		plt.xticks(fontsize=20)
	        plt.ylabel('U-V color',fontsize=40)
	        plt.xlabel('Time (Gyr)',fontsize=40)
	        plt.plot([0,9],[1.3,1.3],'k--')

	plt.legend(loc='upper left',fontsize=22)
	return

def BrightBHGal(bhhalo,bhorbit,filelist='files.list',stepfile='steps.list',dt='100 Myr',lcut=1e43,filename='brightBHgal.pkl'):

	f = open(stepfile,'r')
        steps = np.array(f.readlines()).astype('int')
	f.close()
	f = open(filelist,'r')
	files = f.readlines()
	f.close()

	dtunits = pynbody.units.Unit(dt)
	dt = pynbody.array.SimArray(1.0,dtunits)
	init = np.array([])
	data = {'halonum':init,'step':init,'time':init,'z':init,'sfr':init,'starmass':init,'dmmass':init,'gasmass':init,'metalicity':init,'HIfrac':init,'BHmass':init,'BHlum':init,'interact':init,'offset':init,'otheroffset':init,'otherBHmass':init,'otherBHlum':init,'BHiord':init,'otherBHiord':init,'BHmdot':init,'otherBHmdot':init}
	
	for i in range(len(files)):
                print "gettign data for for ", files[i].strip('\n')
                s = pynbody.load(files[i].strip('\n'))
		if len(bhhalo['mass'][(bhhalo['mass'][:,i]>0)])==0:
			print "no BHs this step! moving on"
                        continue
	#	if len(s.stars[(s.stars['tform']<0)]) == 0: 
#			print "no BHs this step! moving on"
#			continue
                h = s.halos()
		s.physical_units()
		simtime = s.properties['time'].in_units(dtunits)
		
		iords = bhhalo['iord'][(np.in1d(bhhalo['iord'],s.stars['iord'][(s.stars['tform']<0)]))]

		halos = init
		lummean = init
		dist = init
		time = init
		mass = init
		interact = init
		iord = init
		acc = init
		starM = init
		dmM = init
		gasM = init
		
		for id in iords:
			curbh, = np.where(bhorbit['iord']==id)
			curbh2, = np.where(bhhalo['iord']==id)
			curbh = curbh[0]
			curbh2 = curbh2[0]
			if bhhalo['haloID'][curbh2,i]==0: continue
			tt, = np.where((bhorbit['data'][curbh]['Time'].in_units(dtunits) <= simtime)&(bhorbit['data'][curbh]['Time'].in_units(dtunits) >= simtime - dt)&(bhorbit['data'][curbh]['mass'].in_units('Msol')>=1e6))
			lum = bhorbit['data'][curbh]['mdot'][tt].in_units('g s**-1')*0.1*3e10*3e10
			mdot = bhorbit['data'][curbh]['mdot'][tt].in_units('Msol yr**-1')
			tstep = bhorbit['data'][curbh]['dt'][tt].in_units(dtunits)
			lmean = np.sum(lum*tstep)/dt
			mdotmean = np.sum(mdot*tstep)/dt
			if lmean > lcut:
				halos = np.append(halos,bhhalo['haloID'][curbh2,i])
				lummean = np.append(lummean,lmean)
				if bhhalo['distcen'][curbh2,i]>0:
					dist = np.append(dist, bhhalo['distcen'][curbh2,i])
				else:
					dist = np.append(dist, bhhalo['dist'][curbh2,i])
				interact = np.append(interact,bhhalo['interact'][curbh2,i])
				mass = np.append(mass,bhhalo['mass'][curbh2,i])
				iord = np.append(iord,id)
				acc = np.append(acc,mdotmean)
				starM = np.append(starM,bhhalo['halostarmass'][curbh2,i])
				gasM = np.append(gasM,bhhalo['halogasmass'][curbh2,i])
				dmM = np.append(dmM,bhhalo['halodarkmass'][curbh2,i])
		print len(lummean)
		uhalos,invind = np.unique(halos,return_inverse=True)
		data['halonum'] = np.append(data['halonum'],uhalos)
		data['step'] = np.append(data['step'],np.ones(len(uhalos))*steps[i])
		data['time'] = np.append(data['time'],np.ones(len(uhalos))*s.properties['time'].in_units('Gyr'))
		data['z'] = np.append(data['z'],np.ones(len(uhalos))*(s.properties['a']**-1 -1))
		for hh in uhalos:
			print "getting BH info for halo ", hh
			o, = np.where(halos==hh)
			ltmp = lummean[o]
			dtmp = dist[o]
			mtmp = mass[o]
			itmp = iord[o]
			atmp = acc[o]
			smtmp = starM[o]
			gmtmp = gasM[o]
			dmMtmp = dmM[o]
			inttmp = interact[o]
			print ltmp
			lsrt = np.argsort(ltmp)
			data['BHlum'] = np.append(data['BHlum'],ltmp[lsrt[-1]])
			data['offset'] = np.append(data['offset'],dtmp[lsrt[-1]])
			data['BHmass'] = np.append(data['BHmass'],mtmp[lsrt[-1]])
			data['BHiord'] = np.append(data['BHiord'],itmp[lsrt[-1]])
			data['BHmdot'] = np.append(data['BHmdot'],atmp[lsrt[-1]])
			data['starmass'] = np.append(data['starmass'],smtmp[lsrt[-1]])
			data['dmmass'] = np.append(data['dmmass'],dmMtmp[lsrt[-1]])
			data['gasmass'] = np.append(data['gasmass'],gmtmp[lsrt[-1]])
			data['interact'] = np.append(data['interact'],inttmp[lsrt[-1]])

			if len(ltmp) > 1:
				data['otherBHlum'] = np.append(data['otherBHlum'],ltmp[lsrt[-2]])
	                        data['otheroffset'] = np.append(data['otheroffset'],dtmp[lsrt[-2]])
	                        data['otherBHmass'] = np.append(data['otherBHmass'],mtmp[lsrt[-2]])
	                        data['otherBHiord'] = np.append(data['otherBHiord'],itmp[lsrt[-2]])
				data['otherBHmdot'] = np.append(data['otherBHmdot'],atmp[lsrt[-2]])
			else:
				data['otherBHlum'] = np.append(data['otherBHlum'],0)
                                data['otheroffset'] = np.append(data['otheroffset'],0)
                                data['otherBHmass'] = np.append(data['otherBHmass'],0)
                                data['otherBHiord'] = np.append(data['otherBHiord'],0)
				data['otherBHmdot'] = np.append(data['otherBHmdot'],0)

			if len(h[hh].s)>0: 
				newstars, = np.where(h[hh].s['tform'].in_units(dtunits)>=simtime-dt)
				s.read_starlog()
				sfr = h[hh].s['massform'][newstars].in_units('Msol').sum()/dt.in_units('yr')
			else:
				sfr = 0
                        if len(h[hh].g) > 0: 
				metalicity = np.sum(h[hh].g['metals']*h[hh].g['mass'])/np.sum(h[hh].g['mass'])
                        	HI = np.sum(h[hh].g['HI']*h[hh].g['mass'])/np.sum(h[hh].g['mass'])
			else:
				metalicity = 0
				HI = 0
                        data['sfr']=np.append(data['sfr'],sfr)
                        data['metalicity']=np.append(data['metalicity'],metalicity)
                        data['HIfrac'] = np.append(data['HIfrac'],HI)
		del(s)
		del(h)
		gc.collect()
	if filename:
		f = open(filename,'wb')
		pickle.dump(data,f)
		f.close()
	return data
