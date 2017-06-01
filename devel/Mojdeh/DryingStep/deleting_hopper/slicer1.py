import numpy as np 
import math, random, copy, re
from shutil import copyfile 

# The left and right (x) boundaries of the slice
xlo = 20
xhigh = 50

readinFile = "finished_coat.data"
slicedFile = "sliced.data"

copyfile(readinFile,slicedFile)

with open(readinFile,"r") as inFile:
	with open(slicedFile,'w') as outFile:
		#we open the readinfile and the copied sliced file
		foundAtoms = False
		atomCount = 0
		linesToWrite = []
		lines = inFile.readlines() # lines takes every line inFile and turns it into an array
		for line in lines:
			if foundAtoms == False:
				linesToWrite.append(line) #goes to a new line
				if re.match('^[ \t]*Atoms[ \t]*\n',line):
					foundAtoms = True
					linesToWrite.append("\n")		
			elif len(line.split()) > 0: # if it's not a blank line. the spliter consider every cell after apace as a variable in the array
				xpos = float(line.split()[-6])  # we assume that it is the x position
				if xpos >= xlo and xpos <= xhigh:
					atomCount += 1
					# Edit the line
					values = line.split()
					adjustedXPos = xpos - xlo # maintain the original value of x 
					values[0] = str(atomCount)
					values[-6] = str(adjustedXPos)
					currentline = "" #?
					for i in range (len(values)):
						currentline += values[i] + " "
						indeces_to_keep = [0]
					currentline = currentline[:-1] + "\n" # goes to a new line 	
					linesToWrite.append(currentline)
					# currrent line gets redifined every time
				
		# Go back and correct number of atoms
		linesToWrite[1] = str(atomCount) + " atoms\n"

		# Change boundaries of simulation to match the slice size
		linesToWrite[5] = str(0.0000) + " " + str(xhigh - xlo) + " xlo xhi\n"

		for line in linesToWrite:
			outFile.write(line)

	outFile.close()
inFile.close()

print ("writing new file", slicedFile)

