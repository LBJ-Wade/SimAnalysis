import volumeAnalysis
import bhanalysis
import matplotlib.colors as pltcolors
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pynbody
import matplotlib.patches as mpatches
import gc
import readcol

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
				labelstr = False
				if hlabels[hcnt] or simlabels[simcnt]: labelstr = simlabels[simcnt]+" "+hlabels[hcnt]
				plt.scatter(color1[hcnt,:],color2[hcnt,:],c=1./red,norm=redcNorm,cmap=hcmap[hcnt],marker=simmark[simcnt],s=msize,label=labelstr,color=hcol[hcnt])
			hcnt += 1
		if average==True:
			labelstr=False
			if hlabels[0]  or simlabels[0]  : labelstr = simlabels[0]+" "+hlabels[0]
			plt.scatter(color1.mean(axis=0),color2.mean(axis=0),c=1./red,norm=redcNorm,cmap=hcmap[0],marker=simmark[simcnt],s=msize,label=labelstr,color=hcol[0])
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

def colortime(halos=[1,2],average=False,shade=False,simdirs=['romulus8.256gst3.bwBH'],dodust=True,dustfile = 'dust.pkl',magfile = 'mags.pkl',hcol=['red','blue'],simstyle=['solid','solid'],plotData=True,simlabels=['Romulus'],hlabels=['halo 1','halo 2'],overplot=False):
	simcnt = 0
	print "number of halos", len(halos)
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
                                color[hcnt,:] += reddening
			if average==False and shade==False:
				labelstr = None
                                if hlabels[hcnt] or simlabels[simcnt]: labelstr = simlabels[simcnt]+" "+hlabels[hcnt]
				plt.plot(time,color[hcnt,:],color=hcol[hcnt],linewidth=2,linestyle=simstyle[simcnt],label=labelstr)
			hcnt += 1
		if average==True and shade==False:
			labelstr = None
			if hlabels[0] or simlabels[0]: labelstr = simlabels[0]+" "+hlabels[0]
			plt.errorbar(time,color.mean(axis=0),color=hcol[0],yerr=color.std(axis=0),linewidth=2,linestyle=simstyle[simcnt],label=labelstr)
		if shade==True and average==False:
			labelstr = None
                        if hlabels[0] or simlabels[0]: labelstr = simlabels[0]+" "+hlabels[0]
			plt.fill_between(time,color.min(axis=0),color.max(axis=0),facecolor=hcol[0],alpha=0.5)	
		if shade==True and average==True:
			labelstr = None
                        if hlabels[0] or simlabels[0]: labelstr = simlabels[0]+" "+hlabels[0]
			print labelstr
			plt.plot(time,color.mean(axis=0),color=hcol[0],label=labelstr,linestyle=simstyle[simcnt],linewidth=2)
                        plt.fill_between(time,color.mean(axis=0)-color.std(axis=0),color.mean(axis=0)+color.std(axis=0),facecolor=hcol[0],alpha=0.5,label=labelstr)
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

def plotStarGrowth(stars,type='Total',halo=1,color='blue',linestyle='solid',label=None,plotData=True,error=True,datacolor=['grey','k']):
	if plotData==True:
		f = open('/nobackupp8/mtremmel/DATA/Morishita15/Morishita15.pkl','rb')
		data = pickle.load(f)
		f.close()
	if type == 'Total':
		plt.plot(stars['z'],stars['Total'][halo-1,:],color=color,linestyle=linestyle,label=label)
		if plotData==True:
			if error == False: 
				plt.plot(data['z'],data['totMW'],color=datacolor[0],marker='D',label='CANDELS MW',markersize=10,linestyle='')
#				plt.plot(data['z'],data['totMG'],color=datacolor[1],marker='^',label='CANDELS M31',markersize=10,linestyle='')
			if error == True: 
				plt.errorbar(data['z'],data['totMW'],color=datacolor[0],yerr=[data['totMWerrMinus'],data['totMWerrPlus']],fmt='D',label='CANDELS MW',markersize=10)
#				plt.errorbar(data['z'],data['totMG'],color=datacolor[1],yerr=[data['totMGerrMinus'],data['totMGerrPlus']],fmt='^',label='CANDELS M31',markersize=10)
	if type == 'Inner':
                plt.plot(stars['z'],stars['M <2.5 kpc'][halo-1,:],color=color,linestyle=linestyle,label=label)
                if plotData==True:
                        if error == False:
                                plt.plot(data['z'],data['innerMW'],color=datacolor[0],marker='D',label='CANDELS MW',markersize=10,linestyle='')
 #                               plt.plot(data['z'],data['innerMG'],color=datacolor[1],marker='^',label='CANDELS M31',markersize=10,linestyle='')
                        if error == True: 
                                plt.errorbar(data['z'],data['innerMW'],color=datacolor[0],yerr=[data['innerMWerrMinus'],data['innerMWerrPlus']],fmt='D',label='CANDELS MW',markersize=10)
  #                              plt.errorbar(data['z'],data['innerMG'],color=datacolor[1],yerr=[data['innerMGerrMinus'],data['innerMGerrPlus']],fmt='^',label='CANDELS M31',markersize=10)
	if type == 'Outer':
                plt.plot(stars['z'],stars['M >2.5 kpc'][halo-1,:],color=color,linestyle=linestyle,label=label)
                if plotData==True:
                        if error == False:
                                plt.plot(data['z'],data['outerMW'],color=datacolor[0],marker='D',label='CANDELS MW',markersize=10,linestyle='')
   #                             plt.plot(data['z'],data['outerMG'],color=datacolor[1],marker='^',label='CANDELS M31',markersize=10,linestyle='')
                        if error == True: 
                                plt.errorbar(data['z'],data['outerMW'],color=datacolor[0],yerr=[data['outerMWerrMinus'],data['outerMWerrPlus']],fmt='D',label='CANDELS MW',markersize=10)
    #                            plt.errorbar(data['z'],data['outerMG'],color=datacolor[1],yerr=[data['outerMGerrMinus'],data['outerMGerrPlus']],fmt='^',label='CANDELS M31',markersize=10)

	plt.ylabel(r'M$_{*}$ [M$_{\odot}$]',fontsize=30)
	plt.xlabel(r'Redshift',fontsize=30)
	plt.legend(fontsize=20)

	return


def MStarGrowth(lowz,catfiles,Rth='2.5 kpc',halonum=[1],filename='MstarGrowth.pkl'):
        f = open(catfiles,'r')
        files = f.readlines()
        slz = pynbody.load(lowz)
        lz = str(round(slz.properties['a']**-1 -1,3))
	Stars = {'Total':pynbody.array.SimArray(np.zeros((len(halonum),len(files)+1)),'Msol'),
		 'Rhalf':pynbody.array.SimArray(np.zeros((len(halonum),len(files)+1)),'kpc'),
                 'M <'+Rth:pynbody.array.SimArray(np.zeros((len(halonum),len(files)+1)),'Msol'),
		 'M >'+Rth:pynbody.array.SimArray(np.zeros((len(halonum),len(files)+1)),'Msol'),
		 'z':np.zeros(len(files)+1),
		 'halos':halonum}

        catend = '.cat.z'+lz+'\n'
	for i in range(len(files)):
		print "calculating stellar buildup for halos in ", files[i].strip('\n')
		xx = files[i].find(catend)
                simname=files[i][0:xx]
                s = pynbody.load(simname)
                h = s.halos()
		Stars['z'][i] = s.properties['a']**-1 -1
		s.physical_units()
                catf = open(files[i].strip('\n'))
                cat = pickle.load(catf)
                catf.close()
		amigastat = readcol.readcol(simname+'.amiga.stat',asdict=True)
		for j in range(len(halonum)):
                        print "halo", halonum[j]
                        progs, = np.where(cat==halonum[j])
                        if len(progs)==0:
                                print "no progenitors found in this step!"
                                continue
			curhalos, = np.where(np.in1d(amigastat['Grp'],progs))
                        main = amigastat['Grp'][curhalos][np.argmax(amigastat['StarMass(M_sol)'][curhalos])]
                        print "progenitor", main
			h1 = h[main]
			try:
				pynbody.analysis.halo.center(h1,mode='hyb',wrap=True)
			except:
				pynbody.analysis.halo.center(h1,mode='hyb',wrap=True,vel=False)
			Stars['Total'][j,i] = h1.stars['mass'].in_units('Msol').sum()
			Stars['M <'+Rth][j,i] = h1.stars['mass'][(h1.stars['r'].in_units(Rth) <= 1)].in_units('Msol').sum()
			Stars['M >'+Rth][j,i] = h1.stars['mass'][(h1.stars['r'].in_units(Rth) > 1)].in_units('Msol').sum()
			order = np.argsort(h1.stars['r'])
			Mcum = np.cumsum(h1.stars['mass'][order].in_units('Msol'))
			o, = np.where(Mcum >= Stars['Total'][j,i]/2.)
			Stars['Rhalf'][j,i] = h1.stars['r'].in_units('kpc')[order[o[0]]]
			del(Mcum)
			del(order)
			del(o)
			gc.collect()
		del(s)
		del(h)
		del(cat)
		gc.collect()
	h = slz.halos()
        slz.physical_units()
	print "calculating stellar buildup for halos in ",lowz
	Stars['z'][i+1] =  slz.properties['a']**-1 -1
        for j in range(len(halonum)):
                h1 = h[halonum[j]]
		try:
                	pynbody.analysis.halo.center(h1,mode='hyb',wrap=True)
                except:
                        pynbody.analysis.halo.center(h1,mode='hyb',wrap=True,vel=False)
		Stars['Total'][j,i+1] = h1.stars['mass'].in_units('Msol').sum()
                Stars['M <'+Rth][j,i+1] = h1.stars['mass'][(h1.stars['r'].in_units(Rth) <= 1)].in_units('Msol').sum()
                Stars['M >'+Rth][j,i+1] = h1.stars['mass'][(h1.stars['r'].in_units(Rth) > 1)].in_units('Msol').sum()
		order = np.argsort(h1.stars['r'])
                Mcum = np.cumsum(h1.stars['mass'][order].in_units('Msol'))
                o, = np.where(Mcum >= Stars['Total'][j,i+1]/2.)
                Stars['Rhalf'][j,i+1] = h1.stars['r'].in_units('kpc')[order[o[0]]]
	if filename:
		print "saving data..."
		f = open(filename,'wb')
		pickle.dump(Stars,f)
		f.close()
	return Stars

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
			tt, = np.where((bhorbit['data'][curbh]['Time'].in_units(dtunits) <= simtime)&(bhorbit['data'][curbh]['Time'].in_units(dtunits) >= simtime - dt)&(bhorbit['data'][curbh]['mass'].in_units('Msol')-bhorbit['data'][curbh]['dM'].in_units('Msol')>=1e6))
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


def MulaneyBHSFR(sfr):
	#Mulaney+ 2012
	return 1e-3 * sfr

def ChenBHSFR(sfr):
	#Chen+ 2013
	return 10**-3.72 * sfr**1.05

def TotalIRto60micron(LIR):
	'''based on figure 3 of Chary and Elbaz 2001'''
	return LIR/2.

def SFRtoIR(SFR):
	'''based on Daddi+ 2007'''
	return SFR / 1.73e-10 * 3.846e33

def RosarioBHSFR(sfr):
	#Rosario+ 2012
	LIR = SFRtoIR(sfr)
	L60 = TotalIRto60micron(LIR)
	logL60norm = np.log10(L60/1e44)
	logLBH = np.log10(L60/1e44)**(1.0/0.78) + np.log10(6.3e44)
	LBH = 10**logLBH
	mdot = pynbody.array.SimArray(LBH/(0.1*3e10*3e10),'g s**-1')
	return mdot.in_units('Msol yr**-1')

def plotBHSFRdata():
	sfrline = 10**np.arange(-4,4,0.1)
	mdotM12 = MulaneyBHSFR(sfrline)
	mdotC13 = ChenBHSFR(sfrline)
#	mdotR12 = RosarioBHSFR(sfrline)
	plt.plot(sfrline,mdotM12,'b-',label='Mulaney+ 12',linewidth=2)
	plt.plot(sfrline,mdotC13,'b--',label='Chen+ 13',linewidth=2)
#	plt.plot(sfrline,mdotR12,'g-',label='Rosario+ 12',linewidth=2)
	plt.legend(fontsize=20)
	plt.ylabel(r'$\dot{\mathrm{M}_{BH}}$ [M$_{odot}$ yr$^{-1}$]',fontsize=30)
	plt.xlabel(r'SFR [M$_{odot}$ yr$^{-1}$]',fontsize=30)
	return	
