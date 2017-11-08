# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
import random

class DatafileGenerator():
	newFileName = "fcc_packed.data"
	positionLines = []
	velocityLines = []
	xline = []

	# data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
	#cbdStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
	#actStr = "%d %d %d 4.79 0.0 1.0 %f %f %f\n" #string for active particles
	#solStr = "%d %d %d 1.028 0.0 1.0 %f %f %f\n" #string for solvent particles
	#wallStr = "%d %d %d 5.0 0.0 1.0 %f %f %f\n" #string for wall particles
	#cbdStr = "%d %d %d %f %f %f\n" #string for cbd particles
	solStr = "%d %d %f %f %f\n" #string for solvent particles
	#actStr = "%d %d %d %f %f %f\n" #string for active particles
	wallStr = "%d %d %d %f %f %f\n" #string for wall particles

	m_wall = 3.896
   
	m_sol = 4.31      
   
	dia = 2.0

	sol_dia = 2.0
	
	solvent_type,wall_type = 1,2
	rho_sol = 1.028

	solid_d = 3
	xlo = 0.0
	xhi = xlo + int(10.0*dia)
	ylo = 0.0
	yhi = int(10.0*dia)
	zlo = 0.0
	zhi = int(10.0*dia)


	ID = 0
	molID = 0
	def generateDataFile(self):
			self.setUpSimulation()
			self.writeFile()
			print ("generateDataFile is working")

	def setUpSimulation(self):
		xlen = abs(self.xlo-self.xhi)
		ylen = abs(self.yhi-self.ylo)
		zlen = abs(self.zhi-self.zlo)

		dia = self.dia*0.7 # wall diameter
		radius = self.dia/2
		wall_radius = self.dia/3

		vertex1 = (2,0)
		vertex2 = (0.8*xlen,0.7*zlen)

		#Wall far XZ wall
		vertexa = (0,ylen-radius/2,0)
		vertexb = (xlen,ylen,zlen)
		#Wall near XZ wall
		vertexc = (0,0,0)
		vertexd = (xlen,dia,zlen)

		#Wall left YZ wall
		vertexe = (0,wall_radius,0)
		vertexf = (dia,ylen,zlen)

		#Wall right YZ wall
		vertexg = (xlen-dia,wall_radius,0)
		vertexh = (xlen,ylen,zlen)

		#Wall bottom XY wall
		vertexi = (dia,0,0)
		vertexj = (xlen-dia,ylen,dia)
		print ("setUpSimulation is working")

		self.FccVertexes(vertex1,vertex2,self.yhi,1)

		#self.DrawWallVtx(vertexa,vertexb,wall_radius,'XZ')
		self.DrawWallVtx(vertexc,vertexd,wall_radius,'XZ')
		self.DrawWallVtx(vertexe,vertexf,wall_radius,'YZ')
		self.DrawWallVtx(vertexg,vertexh,wall_radius,'YZ')
		#self.DrawWallVtx(vertexi,vertexj,wall_radius,'XY')
		
	def appendLine(self,atomType,x,y,z):
		self.ID += 1
		#if atomType == self.cbd_type:
			#self.positionLines.append(self.cbdStr % (self.ID, self.molID, atomType, x, y, z))
		if atomType == self.solvent_type:
			self.positionLines.append(self.solStr % (self.ID, self.molID, atomType, x, y, z))
		#elif atomType == self.active_type:
			#self.positionLines.append(self.actStr % (self.ID, self.molID, atomType, x, y, z))
		else:
			self.positionLines.append(self.wallStr % (self.ID, self.molID, atomType, x, y, z))

	# FCC packed structure 
	def FccVertexes(self,vertexa,vertexb,yhi,atomType=solvent_type):
			xlo = vertexa[0]
			xhi = vertexb[0]
			zlo = vertexa[1]
			zhi = vertexb[1]
			ylo = self.ylo

			print ("FccVertexes is working")

			self.FCC(xlo,xhi,ylo,yhi,zlo,zhi,atomType)


	def FCC(self,xlo,xhi,ylo,yhi,zlo,zhi,atomType=solvent_type):

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
	def DrawWallVtx(self,Vertex1,Vertex2,radius,plane):
		xlo = Vertex1[0]
		xhi = Vertex2[0]
		ylo = Vertex1[1]
		yhi = Vertex2[1]
		zlo = Vertex1[2]
		zhi = Vertex2[2]

		if plane == 'XY':
			self.DrawWallXY(xlo,xhi,ylo,yhi,zlo,zhi,radius)
		if plane == 'XZ':
			self.DrawWallXZ(xlo,xhi,ylo,yhi,zlo,zhi,radius)
		if plane == 'YZ':
			self.DrawWallYZ(xlo,xhi,ylo,yhi,zlo,zhi,radius)

	def DrawWallXZ(self,xlo,xhi,ylo,yhi,zlo,zhi,radius):
		atomType = self.wall_type
		beta = 0
		dia = 2*radius
		for k in np.arange(zlo,zhi,radius*np.sqrt(3)):
			beta += 1
			for i in np.arange(xlo,xhi,dia):
				x = (-1)**beta*(radius/2) + i
				y = ylo + (yhi-ylo)/2
				z = k
				if x >= 0 and y >= 0 and z >= 0:
					if x <= self.xhi and y <= self.yhi and z <= self.zhi:
						self.appendLine(atomType,x,y,z)

	def DrawWallXY(self,xlo,xhi,ylo,yhi,zlo,zhi,radius):
		atomType = self.wall_type
		beta = 0
		dia = 2*radius
		for k in np.arange(ylo,yhi,radius*np.sqrt(3)):
			beta += 1
			for i in np.arange(xlo,xhi,dia):
				x = (-1)**beta*(radius/2) + i
				y = k
				z = zlo + (zhi-zlo)/2
				if x >= 0 and y >= 0 and z >= 0:
					if x <= self.xhi and y <= self.yhi and z <= self.zhi:
						self.appendLine(atomType,x,y,z)

	def DrawWallYZ(self,xlo,xhi,ylo,yhi,zlo,zhi,radius):
		atomType = self.wall_type
		beta = 0
		dia = 2*radius
		for k in np.arange(zlo,zhi,radius*np.sqrt(3)):
			beta += 1
			for i in np.arange(ylo,yhi,dia):
				y = (-1)**beta*(radius/2) + i
				x = xlo + (xhi-xlo)/2
				z = k
				if x >= 0 and y >= 0 and z >= 0:
					if x <= self.xhi and y <= self.yhi and z <= self.zhi:
						self.appendLine(atomType,x,y,z)	
	def writeFile(self):
		print("writeFile is working")
		print("Writing new file: %s" % self.newFileName)

		xlo = self.xlo
		xhi = self.xhi
		ylo = self.ylo
		yhi = self.yhi
		zlo = self.zlo
		zhi = self.zhi

		with open(self.newFileName,"w") as outFile:

			outFile.write("\n")
			outFile.write("%d atoms\n" % self.ID)
			outFile.write("\n")
			outFile.write("%d atom types\n" % 2)
			outFile.write("\n")
			outFile.write("%f %f xlo xhi\n" % (xlo,xhi))   
			outFile.write("%f %f ylo yhi\n" % (ylo,yhi))
			outFile.write("%f %f zlo zhi\n" % (zlo,zhi))

			outFile.write("\n")
			outFile.write("Masses\n")
			outFile.write("\n")
			outFile.write("1 %f\n" % self.m_sol)
			outFile.write("2 %f\n" % self.m_wall)

			outFile.write("\n")
			outFile.write("Atoms\n")
			outFile.write("\n")

			for line in self.positionLines:
				outFile.write(line)
			for line in self.xline:
				outFile.write(line)

if __name__ == "__main__":
	generator = DatafileGenerator()
	generator.generateDataFile()

