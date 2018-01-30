# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
#from PIL import Image
import random, math, copy

class DatafileGenerator():
	newFileName = "fcc+wall.data"
	positionLines = []

   # data string containing atom ID, molecule ID, type, rho, e, cv, x, y, z 
	#cbdStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
	#solStr = "%d %d %d 1.028 0.0 1.0 %f %f %f\n" #string for solvent particles
	#actStr = "%d %d %d 4.79 0.0 1.0 %f %f %f\n" #string for active particles
	#wallStr = "%d %d %d 2 0.0 1.0 %f %f %f\n" #string for wall particles
	cbdStr = "%d %d %d %f %f %f\n" #string for cbd particles
	solStr = "%d %d %d %f %f %f\n" #string for solvent particles
	actStr = "%d %d %d %f %f %f\n" #string for active particles
	wallStr = "%d %d %d %f %f %f\n" #string for wall particles

	m_cbd = 3.896
	m_sol = 5.36
	m_act = 20.064
	m_wall = 1.047

	dia = 2.0
	cbd_dia = 2.0
	act_dia = 8.0
	sol_dia = 2.0
	wall_dia = 2.0

	cbd_type,solvent_type,active_type,wall_type = 1,2,3,4

	active_thickness_factor = 1

	True_activecount = 0
	activecount = 0
	solventcount = 0
	cbdcount = 0

	'''activeShapes = ['sphere','ellipsoid','rod','single','diparticle']
	activeShape = activeShapes[4]
	if activeShape == 'sphere':
		active_dia = 8.0
	elif activeShape == 'ellipsoid':
		active_dia = 18.0
	elif activeShape == 'rod':
		active_dia = 12.0
	elif activeShape == 'single':
		active_dia = 3.0
	elif activeShape == 'diparticle': 
		active_dia = 4.0'''

	scale = 1

	shrink = 1

	xlo = 0.0
	xhi = int(100.0*dia)
	ylo = 0.0
	yhi = int(100.0*dia)
	zlo = 0.0
	zhi = int(300.0*dia)

	ID = 0
	molID = 1

	def generateDataFile(self):
		self.setUpSimulation()
		self.writeFile()

	def setUpSimulation(self):
		xlen = abs(self.xhi-self.xlo)
		ylen = abs(self.yhi-self.ylo)
		zlen = abs((self.zhi-self.zlo)/self.shrink)

		vertexalpha = (.4*xlen,zlen*.05 + self.dia)
		vertexbeta =  (.6*xlen,zlen*.9)

		vertex9 = (0,0,zlen*0.05/2)
		vertex10 = (xlen,ylen,zlen*0.05/2)

		vertex1 = (2,1.5)
		vertex2 = (0.8*xlen,zlen*0.7)
		vertex4 = (xlen,zlen*0.97)

		
	  
		#self.drawWallFromVtxs(vertex9,vertex10)

		#self.FccVertexes(vertex1,vertex2,self.yhi,1)
		
		self.fillCubeWithRandomMixVtxs(vertex1,vertex4)
		#self.fillCubeWithCBDVtxs(vertex5,vertex6,checkForOverlap=False)



	'''def fillCubeWithRandomMixVtxs(self,v1,v2):
		(x1,z1) = v1
		(x2,z2) = v2
		self.fillCubeWithRandomMix(self.xlo,self.ylo,self.zhi,x2,self.yhi,z2)

	def fillCubeWithRandomMix(self,x,y,z,x2,y2,z2):
		xl = min(x,x2)
		xh = max(x,x2)
		yl = min(y,y2)
		yh = max(y,y2)
		zl = min(z,z2)
		zh = max(z,z2)
		radius = self.act_dia/2
		for yi in np.arange(yl+self.dia*3,yh-self.dia,self.dia*4):
			for zi in np.arange(zl--self.dia*5,zh-self.dia*7,self.dia*4):
				for xi in np.arange(xl+self.dia*3,xh-self.dia*2,self.dia*4):
					val = random.randint(0,150)
					if val == 0:
						self.activecount += 1
						self.drawRaspberry(xi,yi,zi,radius)
						
					elif val >= 1 and val < 40:
						self.cbdcount += 1
						atom_type = self.cbd_type
						#self.FCC1(xlo,xhi,ylo,yhi,zlo,zhi,atomType=cbd_type)
						self.appendLine(atom_type,xi+self.dia*5/2,yi+self.dia/2,zi+self.dia/2)
					elif val >= 40 and val < 151:
						self.solventcount += 1
						atom_type = self.solvent_type
						#self.FCC2(xlo,xhi,ylo,yhi,zlo,zhi,atomType=solvent_type)
						self.appendLine(atom_type,xi+self.dia*1/2,yi+self.dia/2,zi+self.dia/2)'''   
	

	def FccVertexes(self,vertexa,vertexb,yhi,atomType):
			xlo = vertexa[0]
			xhi = vertexb[0]
			zlo = vertexa[1]
			zhi = vertexb[1]
			ylo = self.ylo

			print ("FccVertexes is working")

			self.FCC(xlo,xhi,ylo,yhi,zlo,zhi,atomType)


	def FCC1(self,xlo,xhi,ylo,yhi,zlo,zhi,atomType=cbd_type):

		for k in np.arange(zlo-zhi,5*zhi,self.dia/2):
			for j in np.arange(ylo-yhi,5*yhi,self.dia/2):
				for i in np.arange(xlo-xhi,5*xhi,self.dia/2):
					x = (2*i + (( j + k )%2))
					y = (np.sqrt(3)*( j + (1/3)*(k%2)))
					z = ((2*np.sqrt(6)/3)*k)
					if z > zhi or x > xhi or y > yhi:
						break
					if z > zlo:
						if x > xlo:
							if y > ylo:
								self.appendLine(atomType,x,y,z)
	
	def FCC2(self,xlo,xhi,ylo,yhi,zlo,zhi,atomType=solvent_type):

		for k in np.arange(zlo-zhi,5*zhi,self.dia/2):
			for j in np.arange(ylo-yhi,5*yhi,self.dia/2):
				for i in np.arange(xlo-xhi,5*xhi,self.dia/2):
					x = (2*i + (( j + k )%2))
					y = (np.sqrt(3)*( j + (1/3)*(k%2)))
					z = ((2*np.sqrt(6)/3)*k)
					if z > zhi or x > xhi or y > yhi:
						break
					if z > zlo:
						if x > xlo:
							if y > ylo:
								self.appendLine(atomType,x,y,z)
	     

	def drawRaspberry(self,x,y,z,radius):
		self.molID += 1
		for xi in np.arange(x-radius,x+radius,self.dia):
			for yi in np.arange(y-radius,y+radius,self.dia):
				for zi in np.arange(z-radius,z+radius,self.dia):
					if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**1.8:
						#if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) > radius**1:  # making hollow
							self.appendLine(self.active_type,xi,yi,zi)

	def appendLine(self,atomType,x,y,z):
		self.ID += 1
		if atomType == self.cbd_type:
			self.positionLines.append(self.cbdStr % (self.ID, self.molID, atomType, x, y, z))
		elif atomType == self.solvent_type:
			self.positionLines.append(self.solStr % (self.ID, self.molID, atomType, x, y, z))
		elif atomType == self.active_type:
			self.positionLines.append(self.actStr % (self.ID, self.molID, atomType, x, y, z))
		else:
			self.positionLines.append(self.wallStr % (self.ID, self.molID, atomType, x, y, z))

	def drawWallFromVtxs(self,vtx1,vtx2,atomType=wall_type):
		(x,y,z)   = vtx1
		(x2,y2,z2) = vtx2
		self.drawWall(x,z,x2,z2,atomType)

	def drawWall(self,x,z,x2,z2,atomType=wall_type):
		self.molID += 1
		print("molID for wall is:", self.molID)
		length = math.sqrt((x2-x)**2 + (z2-z)**2)
		numPoints = math.ceil(length/(self.dia*1))
		delta_x = (x2-x)/numPoints
		delta_z = (z2-z)/numPoints
		for i in range(numPoints):
			xi = x+i*delta_x
			zi = z+i*delta_z
			for yi in np.arange(self.ylo,self.yhi,self.dia*0.25):
				self.appendLine(atomType,xi,yi,zi)

	def drawWallFromVtxs2(self,vtx1,vtx2,atomType=wall_type):
		(x,y,z)   = vtx1
		(x2,y2,z2) = vtx2
		self.drawWall2(x,y,z,x2,y2,z2,atomType)

	def drawWall2(self,x,y,z,x2,y2,z2,atomType=wall_type):
		self.molID += 1
		print("molID for wall is:", self.molID)
		length = math.sqrt((x2-x)**2 + (z2-z)**2)
		numPoints = math.ceil(length/(self.dia*.5))
		delta_x = (x2-x)/numPoints
		delta_y = (y2-y)/numPoints
		delta_z = (z2-z)/numPoints
		for i in range(numPoints):
			xi = x+i*delta_x
			zi = z+i*delta_z
			for yi in np.arange(self.y0,self.y1,self.dia*1):
				self.appendLine(atomType,xi,yi,zi)

  
  
	def writeFile(self):
		print("Writing new file: %s" % self.newFileName)

		xmin = self.xlo
		xmax = self.xhi
		ymin = self.ylo
		ymax = self.yhi
		zmin = self.zlo
		zmax = self.zhi

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
			outFile.write("2 %f\n" % self.m_sol)
			outFile.write("3 %f\n" % self.m_act)
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