# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
#from PIL import Image
import random, math, copy

class DatafileGenerator():
	newFileName = "slurry.data"
	positionLines = []

   # data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
	cbdStr = "%d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
	actStr = "%d %d 2 0.0 1.0 %f %f %f\n" #string for active particles
	solStr = "%d %d 1.028 0.0 1.0 %f %f %f\n" #string for solvent particles
	wallStr = "%d %d 2.0 0.0 1.0 %f %f %f\n" #string for wall particles
	#cbdStr = "%d %d %d %f %f %f\n" #string for cbd particles
	#solStr = "%d %d %d %f %f %f\n" #string for solvent particles
	#actStr = "%d %d %d %f %f %f\n" #string for active particles
	#wallStr = "%d %d %d %f %f %f\n" #string for wall particles

	m_cbd = 3.896
	m_sol = 4.31
	m_act = 8.3 #20.064 #me before:     1284.11       
	m_wall = 28.378    #me before: 2.61

	dia = 2.0
	radius = 1.0

	cbd_dia = 2.0
	act_dia = 2 #????is it correct?
	sol_dia = 2.0
	wall_dia = 3.0

	cbd_type,active_type,solvent_type,wall_type = 1,2,3,4

	active_thickness_factor = 1

	True_activecount = 0
	activecount = 0
	solventcount = 0
	cbdcount = 0


	scale = 0.5

	shrink = 1

	x0 = 0.0
	x1 = int(20.0*dia)
	y0 = 0.0
	y1 = int(10.0*dia)
	z0 = 0
	z1 = int(23.0*dia)
	#z1 = scale*4*50.0*dia*shrink


	ID = 0
	#molID = 0

	def generateDataFile(self):
		self.setUpSimulation()
		self.writeFile()

	def setUpSimulation(self):
		xlen = abs(self.x1-self.x0)
		ylen = abs(self.y1-self.y0)
		zlen = abs((self.z1-self.z0)/self.shrink)

	

		vertex1 = (0.0*xlen,ylen*0.0,zlen*0.0)
		vertex2 = (0.8*xlen,zlen*0.7)
		vertex4 = (1.0*xlen,ylen*1.0,zlen*1.0)
		

		self.fillCubeWithRandomMixVtxs(vertex1,vertex4)
	   

	def fillCubeWithRandomMixVtxs(self,v1,v2):
		(x1,y1,z1) = v1
		(x2,y2,z2) = v2
		self.fillCubeWithRandomMix(x1,y1,z1,x2,y2,z2)
	   

		 

	def fillCubeWithRandomMix(self,x,y,z,x2,y2,z2):
		xl = min(x,x2)
		xh = max(x,x2)
		yl = min(y,y2)
		yh = max(y,y2)
		zl = min(z,z2)
		zh = max(z,z2)
		radius = self.act_dia*2
		for k in np.arange(zl-zh,5*zh,self.dia*0.5):
			for j in np.arange(yl-yh,5*yh,self.dia*0.5):
				for i in np.arange(xl-xh,5*xh,self.dia*0.5):
					xi = (2*i + (( j + k )%2))
					yi = (np.sqrt(3)*( j + (1/3)*(k%2)))
					zi = ((2*np.sqrt(6)/3)*k)
					if zi > zh or xi > xh or yi > yh:
						break
					if zi > zl:
						if xi > xl:
							if yi > yl:
								val = random.randint(0,150)
								if val == 0:
									self.activecount += 1
									self.drawRaspberry(xi,yi,zi,radius)
									
								elif val >= 1 and val < 40:
									self.cbdcount += 1
									atom_type = self.cbd_type
									#self.FccVertexes(self,vertex1,vertex4,yhi,atomType)
									self.appendLine(atom_type,xi+self.dia*5/2,yi+self.dia/2,zi+self.dia/2)
								elif val >= 40 and val < 151:
									self.solventcount += 1
									atom_type = self.solvent_type
									self.appendLine(atom_type,xi+self.dia*1/2,yi+self.dia/2,zi+self.dia/2)
															


	def appendLine(self,atomType,x,y,z):
		self.ID += 1
		if atomType == self.cbd_type:
			self.positionLines.append(self.cbdStr % (self.ID, atomType, x, y, z))
		elif atomType == self.solvent_type:
			self.positionLines.append(self.solStr % (self.ID, atomType, x, y, z))
		elif atomType == self.active_type:
			self.positionLines.append(self.actStr % (self.ID, atomType, x, y, z))
		else:
			self.positionLines.append(self.wallStr % (self.ID, atomType, x, y, z))

	

	def drawRaspberry(self,x,y,z,radius):
		
		for xi in np.arange(x-radius,x+radius,self.dia):
			for yi in np.arange(y-radius,y+radius,self.dia):
				for zi in np.arange(z-radius,z+radius,self.dia):
					if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**1.6:
						#if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) > radius**1:  # making hollow
							self.appendLine(self.active_type,xi,yi,zi)




	def writeFile(self):
		print("Writing new file: %s" % self.newFileName)

		xmin = self.x0
		xmax = self.x1
		ymin = self.y0
		ymax = self.y1
		zmin = self.z0
		zmax = self.z1

		with open(self.newFileName,"w") as outFile:

			outFile.write("\n")
			outFile.write("%d atoms\n" % self.ID)
			outFile.write("\n")
			outFile.write("%d atom types\n" % 8)
			outFile.write("\n")
			outFile.write("%f %f xlo xhi\n" % (xmin,xmax))   
			outFile.write("%f %f ylo yhi\n" % (ymin,ymax))
			outFile.write("%f %f zlo zhi\n" % (zmin,zmax))

			outFile.write("\n")
			outFile.write("Masses\n")
			outFile.write("\n")
			outFile.write("1 %f\n" % self.m_cbd)
			outFile.write("2 %f\n" % self.m_act)
			outFile.write("3 %f\n" % self.m_sol)
			outFile.write("4 %f\n" % self.m_wall)
			outFile.write("5 %f\n" % self.m_wall)
			outFile.write("6 %f\n" % self.m_wall)
			outFile.write("7 %f\n" % self.m_wall)
			outFile.write("8 %f\n" % self.m_wall)

			outFile.write("\n")
			outFile.write("Atoms\n")
			outFile.write("\n")

			for line in self.positionLines:
				outFile.write(line)
	
if __name__ == "__main__":
	generator = DatafileGenerator()
	generator.generateDataFile()