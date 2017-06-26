import numpy as np 
import math, random, copy, re
from shutil import copyfile 

# The left and right (x) boundaries of the slice
xlo = 10
xhigh = 20

readinFile = "finished_coat.data"
slicedFile = "slicedNew.data"

copyfile(readinFile,slicedFile)

with open(readinFile,"r") as inFile:
	with open(slicedFile,'w') as outFile:
		
		foundAtoms = False
		atomCount = 0
		linesToWrite = []
		lines = inFile.readlines()
		for line in lines:
			if foundAtoms == False:
				linesToWrite.append(line)
				if re.match('^[ \t]*Atoms[ \t]*\n',line):
					foundAtoms = True
					linesToWrite.append("\n")		
			elif len(line.split()) > 0:
				xpos = float(line.split()[-6])  # we assume that it is the x position
				if xpos >= xlo and xpos <= xhigh:
					atomCount += 1
					# Edit the line
					values = line.split()
					adjustedXPos = xpos - xlo
					values[0] = str(atomCount)
					values[-6] = str(adjustedXPos)
					currentline = ""
					for i in range (len(values)):
						currentline += values[i] + " "
						indeces_to_keep = [0]
					currentline = currentline[:-1] + "\n"	
					linesToWrite.append(currentline)
					
				
		# Go back and correct number of atoms
		linesToWrite[1] = str(atomCount) + " atoms\n"

		# Change boundaries of simulation to match the slice size
		linesToWrite[5] = str(0.0000) + " " + str(xhigh - xlo) + " xlo xhi\n"

		for line in linesToWrite:
			outFile.write(line)

	outFile.close()
inFile.close()

print ("writing new file", slicedFile)

