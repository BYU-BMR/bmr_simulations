import numpy as np 
import math, random, copy, re
from shutil import copyfile 

#data string containing molecule ID, type, walldia, rho, x, y, z, 0 0 0
wallStr = "%d %d %d 0.93 0.0 1.0 %f %f %f/n" #string for cbd particles


walldia = 2
ID = 0
molID = 0
wall_type = 4

def appendLine(wallid,atomType,x,y,z):
	linesToWrite.append(wallStr % (wallid,molID,atomType,x,y,z))


readinFile = "sliced.data"
preparedFile = "ready_to_calendar.data"

copyfile(readinFile,preparedFile)

with open(readinFile,"r") as inFile:
	with open(preparedFile,'w') as outFile:

		foundAtoms = False
		buildwall = False
		copyatoms = True
		atomCount = 0
		linesToWrite = []
		lines = inFile.readlines()
		

		for line in lines:

			if foundAtoms == False:
				linesToWrite.append(line)
				if re.match('^[ \t]*Atoms.*\n',line):
					foundAtoms = True
					buildwall = True
					linesToWrite.append("\n")

			elif buildwall == True:
				xline = linesToWrite[5]
				yline = linesToWrite[6]
				zline = linesToWrite[7]
				xhi = float(xline.split()[1])
				yhi = float(yline.split()[1])
				ylo = float(yline.split()[0])
				plateheight = float(zline.split()[1]) - walldia
				platelen = float(xhi)
				vertex1 = (0,plateheight)
				vertex2 = (platelen,plateheight)
				buildwall = True
				
				(x,z)   = vertex1
				(x2,z2) = vertex2
				atomType = wall_type
				molID += 1
				print("molID for wall is:", molID)
				numpoints = math.ceil(platelen/walldia)
				delta_x = (x2-x)/numpoints
				for i in range(numpoints):
					xi = x + i*delta_x
					zi = z
					y1 = float(yhi) - walldia/2
					for yi in np.arange(ylo,yhi,walldia):
						ID += 1
						appendLine(ID,atomType,xi,yi,zi)
						buildwall = False
			elif copyatoms == True:
				atomCount = ID
				if re.match('^[ \t]*Atoms.*\n\n',line):
					atomCount += 1
					values = line.split()
					values[0] = str(atomCount)
					currentline = ""
					for value in values:
						currentline += value + " "
					currentline = currentline[:-1] + "\n"
					linesToWrite.append(currentline)
					if re.match('^[ \t]*Velocities.*\n',line):
						copyatoms = False

		# Go back and correct number of atoms
		linesToWrite[1] = str(atomCount) + " atoms\n"

		for line in linesToWrite:
			outFile.write(line)

	outFile.close()
inFile.close()

print ("writing new file", preparedFile)


						
					







				



		 