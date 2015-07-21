#!/usr/bin/python

# This script reads the output files and gets BH info
# this happens without repeating steps

import sys

files=[]
for lines in open('outfiles.list').readlines(): 
	fields=lines.split()
	files.append(fields[0])
	# fields[0] is the first element in the one-element list.

outfile=open('out.bh','w')  # file for outputs

for i in range(len(files)):
	print i, 'file = ', files[i]
	if i != len(files)-1: 
		nextfile=str(files[i+1]) # I want it to read the next file
		# now find the line to stop at.
		f = open(nextfile)
	        for line in f:
        	        if line.startswith('Restart Step:'): 
				N = line[13:-1]+'.00000' # dont include the \n
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






          




	
