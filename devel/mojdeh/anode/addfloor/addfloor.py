import numpy as np 
import math, random, copy, re
from shutil import copyfile 

#data string containing molecule ID, type, walldia, rho, x, y, z, 0 0 0
wallStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
#allStr = "%d %d %d %f %f %f\n"

atomsline = 2
m_wall = 3.72
walldia = 2
ID = 0
wallmolID = 200
wall_type = 5

def appendLine(wallid,atomType,x,y,z):
	linesToWrite.append(wallStr % (wallid,molID,atomType,x,y,z))

readinFile = "finished_disc_shrink_33.data"
preparedFile = "readry.data"

copyfile(readinFile,preparedFile)

with open(readinFile,"r") as inFile:
	with open(preparedFile,'w') as outFile:

		foundAtoms = False
		buildwall = False
		linesToWrite = []
		lines = inFile.readlines()
		
		for line in lines:
			linesToWrite.append(line)
			buildwall = True
			if re.match('^Velocities.*\n',line):
				linesToWrite.pop(-1)
				linesToWrite.pop(-1)
				break

		atoms = linesToWrite[atomsline]
		atomCount = float(atoms.split()[0])
		ID = 10000
		if buildwall == True:
			xline = linesToWrite[5]
			yline = linesToWrite[6]
			zline = linesToWrite[7]
			
			xhi = float(xline.split()[1])
			xlo = float(xline.split()[0])
			yhi = float(yline.split()[1])
			ylo = float(yline.split()[0])
			zlo = float(zline.split()[0]) - 3*walldia
			zhi = float(zline.split()[1])
			
			plateheight = zlo + 3*walldia
			platelen = float(xhi) - float(xlo)
			
			vertex1 = (0,plateheight)
			vertex2 = (platelen,plateheight)
			
			(x,z)   = vertex1
			(x2,z2) = vertex2
			atomType = wall_type
			molID = wallmolID
			print("molID for wall is:", molID)
			numpoints = math.ceil(platelen/walldia)
			delta_x = (x2-x)/numpoints
			
			for i in range(numpoints):
				xi = x + i*delta_x
				zi = z
				y1 = float(yhi) - walldia/2
				for yi in np.arange(ylo,yhi,walldia/3):
					ID += 1
					atomCount += 1
					appendLine(ID,atomType,xi,yi,zi)
					buildwall = False

		# Go back and correct number of atoms
		linesToWrite[atomsline] = str(int(atomCount)) + " atoms\n"
		linesToWrite[7] = str(zlo) + " " + str(zhi) + " zlo" + " zhi\n"
		linesToWrite[14] = "4 "+ str(m_wall) + "\n"
		
		for line in linesToWrite:
			outFile.write(line)

	outFile.close()
inFile.close()

print ("writing new file", preparedFile)


						
					







				



		 