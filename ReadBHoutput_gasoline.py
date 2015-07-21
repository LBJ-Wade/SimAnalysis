#!/usr/bin/python

# This script reads the output files and gets BH info
# this happens without repeating steps

import sys
import numpy as np

files=[]
for lines in open('outfiles.list').readlines(): 
	fields=lines.split()
	files.append(fields[0])
	# fields[0] is the first element in the one-element list.

outfile=open('out.bh','w')  # file for outputs
out_dm = np.zeros((1,5),dtype=('S20'))
out_mdot = np.zeros((1,5),dtype=('S20'))
start=0
start2=0
for i in range(len(files)):
	print i, 'file = ', files[i]
	if i != len(files)-1: 
		nextfile=str(files[i+1]) # I want it to read the next file
		# now find the line to stop at.
		f = open(nextfile)
	        for line in f:
        	        if line.startswith('Restart Step:'): 
				N = line[13:-1]+'.00000'
                        	break
			
	else : N = 'gobbledygook'        

	badline='Calculating Gravity, Step:'+N
	filestring=files[i]+'\n'  
	outfile.write(filestring) # want to print which outfile the lines are from
	f = open(files[i])
	for line in f:  
		if line.startswith(badline): break  #if this line is encountered, quit!
 		else:
			if line.startswith('BHSink'):
					#sweet now I just need to get it in the file.
					outfile.write(line)
					tmp = np.array(line.split())
					if np.size(np.where(tmp=='dm:')) and np.size(np.where(tmp=='dE')):
						tmp[1] = tmp[1].strip(':')
						if start<1: 
							out_dm[0] = tmp[[1,3,5,7,9]]
							start = 1
						if start==1: out_dm = np.append(out_dm,[tmp[[1,3,5,7,9]]],axis=0)
					if np.size(np.where(tmp=='mdot')):
						tmp[1] = tmp[1].strip(':')
						if start2<1: 
							out_mdot[0] = tmp[[1,3,6,9,11]]
							start2=1
						if start2==1: out_mdot = np.append(out_mdot,[tmp[[1,3,6,9,11]]],axis=0)
print out_dm
print out_mdot
out_mdot = out_mdot.astype(np.float)
out_dm = out_dm.astype(np.float)
ff1 = open("out.dm.npy",'w')
np.save(ff1,out_dm)
ff2 = open("out.mdot.npy",'w')
np.save(ff2,out_mdot)
ff1.close()
ff2.close()



	
