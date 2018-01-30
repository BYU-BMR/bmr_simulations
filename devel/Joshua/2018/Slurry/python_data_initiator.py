import numpy as np 
import math, random, copy, re
from shutil import copyfile 

class DatafileGenerator():
	newFileName = "slurry_FCC.data"
	positionLines = []

	# atom data string containing atom_id atom_type rho ev cp x y z
	solventstr = "%d %d 1.43 0.0 1.0 %f %f %f\n"
	#solventstr = "%d %d %f %f %f\n" # use this string instead if previewing in ovito
	cbdstr = "%d %d 2.00 0.0 1.0 %f %f %f\n"
	#cbdstr = "%d %d %f %f %f\n" # use this string instead if previewing in ovito
	wallstr = "%d %d 2.00 0.0 1.0 %f %f %f\n"
	#wallstr = "%d %d %f %f %f\n" # use this string instead if previewing in ovito

	solvent_type,cbd_type,wall_type = 1,2,3

	FCC_latice = True

	dia = 2.0 # micrometer
	density_of_solvent = 1.03 # picogram/micrometer**3
	density_of_cbd = 2.00 # picogram/micrometer**3
	volume = (1/6)*np.pi*dia**3

	m_solvent = density_of_solvent*volume
	m_cbd = density_of_cbd*volume
	m_wall = 2*m_cbd

	print ("mass of solvent =", m_solvent,"picogram/micrometer^3")

	ID = 0

	xlo = 0
	xhi = 20.010000000000003e+00
	ylo = 0
	yhi = 10.009999999999994e+00
	zlo = 0
	zhi = 40.010000000000003e+00


	def generateDataFile(self):
		self.setUpSimulation()
		self.writeFile()
		print ("generateDataFile is working")

	def setUpSimulation(self):
		if self.FCC_latice == True:
			Grid_Gap = 0
		else:
			Grid_Gap = self.dia/2
		xlen = abs(self.xlo-self.xhi)
		ylen = abs(self.yhi-self.ylo)
		zlen = abs(self.zhi-self.zlo)

		dia = self.dia*0.7 # wall diameter
		radius = self.dia/2
		wall_radius = dia/4

		# Water Column
		vertex1 = (0.25*xlen,dia+.2*dia+Grid_Gap)
		vertex2 = (0.75*xlen,0.3*zlen)

		#Wall far XZ wall
		vertexa = (0,ylen-dia,0)
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
		vertexi = (0,0,0)
		vertexj = (xlen,ylen,.75*dia)

		#Position of Wedge
		vertexw = (0.5*xlen - 10*dia,0,0)
		vertexu = (0.5*xlen + 10*dia,0,10*dia)

		print ("setUpSimulation is working")

		self.FccVertexes(vertex1,vertex2,0.65*self.yhi,1)
		#self.DrawWallVtx(vertexa,vertexb,wall_radius,'XZ')
		#self.DrawWallVtx(vertexc,vertexd,wall_radius,'XZ')
		#self.DrawWallVtx(vertexe,vertexf,wall_radius,'YZ')
		#self.DrawWallVtx(vertexg,vertexh,wall_radius,'YZ')
		#self.DrawWallVtx(vertexi,vertexj,wall_radius,'XY')
		#self.DrawWedgeVtx(vertexw,vertexu,wall_radius)

	def appendLine(self,atomType,x,y,z):
		self.ID += 1
		if atomType == self.solvent_type:
			self.positionLines.append(self.solventstr % (self.ID, atomType, x, y, z))
		elif atomType == self.cbd_type:
			self.positionLines.append(self.cbdstr % (self.ID, atomType, x, y, z))
		else:
			self.positionLines.append(self.wallstr % (self.ID, atomType, x, y, z))

	def drawWreckingBall(self,x,y,z,radius):
		self.molID += 1
		for xi in np.arange(x-radius,x+radius,self.dia):
			for yi in np.arange(y-radius,y+radius,self.dia):
				for zi in np.arange(z-radius,z+radius,self.dia):
					if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**2:
						self.appendLine(self.wall,xi,yi,zi)


	def FccVertexes(self,vertexa,vertexb,yhi,atomType):

		xlo = vertexa[0]
		xhi = vertexb[0]
		zlo = vertexa[1]
		zhi = vertexb[1]
		ylo = .25*self.yhi

		print ("FccVertexes is working")

		if self.FCC_latice == True:
			self.FCC(xlo,xhi,ylo,yhi,zlo,zhi,atomType)
		else:
			self.Grid(xlo,xhi,ylo,yhi,zlo,zhi,atomType)


	def FCC(self,xlo,xhi,ylo,yhi,zlo,zhi,atomType):

		for k in np.arange(zlo-2*zhi,2*zhi,self.dia/2):
			for j in np.arange(ylo-2*yhi,2*yhi,self.dia/2):
				for i in np.arange(xlo-2*xhi,2*xhi,self.dia/2):
					x = (2*i + (( j + k )%2))
					y = (np.sqrt(3)*( j + (1/3)*(k%2)))
					z = ((2*np.sqrt(6)/3)*k)
					if z > zhi or x > xhi or y > yhi:
						break
					if z >= zlo:
						if x >= xlo:
							if y >= ylo:
								val = random.randrange(0,100)
								if val >= 0 and val < 50:
									atomType = self.cbd_type
									self.appendLine(atomType,x,y,z)
								else:
									atomType = self.solvent_type
									self.appendLine(atomType,x,y,z)
		print ("FCC is working")

	def Grid(self,xlo,xhi,ylo,yhi,zlo,zhi,atomType):

		for k in np.arange(zlo,zhi,self.dia):
			for j in np.arange(ylo,yhi,self.dia):
				for i in np.arange(xlo,xhi,self.dia):
					x = i
					y = j
					z = k
					val = random.randrange(0,100)
					if val >= 0 and val < 50:
						atomType = self.cbd_type
						self.appendLine(atomType,x,y,z)
					else:
						atomType = self.solvent_type
						self.appendLine(atomType,x,y,z)

		print ("Grid is working")

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
		alpha = -1
		beta = 0
		dia = 2*radius
		for j in np.arange(zlo,zhi,radius*np.sqrt(3)):
			alpha += 1
			for k in np.arange(ylo,yhi,radius*np.sqrt(3)):
				beta += 1
				for i in np.arange(xlo,xhi,dia/2):
					x = (-1)**beta*(radius/4) + i + (-1)**alpha*(radius/8)
					y = k
					z = j
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

	def DrawWedgeVtx(self,Vertex1,Vertex2,radius):
		xlo = Vertex1[0]
		xhi = Vertex2[0]
		ylo = Vertex1[1]
		yhi = Vertex2[1]
		zlo = Vertex1[2]
		zhi = Vertex2[2]

		self.DrawWedge(xlo,xhi,ylo,yhi,zlo,zhi,radius)

	def DrawWedge(self,xlo,xhi,ylo,yhi,zlo,zhi,radius):
		atomType = self.wall_type
		alpha = -1
		beta = 0
		dia = 2*radius
		for j in np.arange(zlo,zhi,radius*np.sqrt(3)):
			alpha += 1
			for k in np.arange(ylo,yhi,radius*np.sqrt(3)):
				beta += 1
				for i in np.arange(xlo,xlo + (xhi-xlo)/2,radius*np.sqrt(3)):
					x = i
					y = k
					z = j
					a = x - xlo
					if a == z:
						if y >= 0:
							if y <= self.yhi:
								self.appendLine(atomType,x,y,z)
		for j in np.arange(zlo,zhi,radius*np.sqrt(3)):
			alpha += 1
			for k in np.arange(ylo,yhi,radius*np.sqrt(3)):
				beta += 1
				for i in np.arange(xlo + (xhi-xlo)/2,xhi,radius*np.sqrt(3)):
					x = i
					y = k
					z = j
					a = x - xlo
					b = z*(-1)
					if a == b:
						if y >= 0:
							if y <= self.yhi:
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
			outFile.write("%d atom types\n" % 3)
			outFile.write("\n")
			outFile.write("%f %f xlo xhi\n" % (xlo,xhi))   
			outFile.write("%f %f ylo yhi\n" % (ylo,yhi))
			outFile.write("%f %f zlo zhi\n" % (zlo,zhi))

			outFile.write("\n")
			outFile.write("Masses\n")
			outFile.write("\n")
			outFile.write("1 %f\n" % self.m_solvent)
			outFile.write("2 %f\n" % self.m_cbd)
			outFile.write("3 %f\n" % self.m_wall)


			outFile.write("\n")
			outFile.write("Atoms\n")
			outFile.write("\n")

			for line in self.positionLines:
				outFile.write(line)

if __name__ == "__main__":
	generator = DatafileGenerator()
	generator.generateDataFile()
